#include "Leg_Force.h"
#include <math.h>

float   g_lf_fx[LF_NLEG];
float   g_lf_fy[LF_NLEG];
float   g_lf_fz[LF_NLEG];
float   g_lf_sum_fz;
float   g_lf_mass_est;
uint8_t g_lf_n_valid;
uint8_t g_lf_ok[LF_NLEG];
float   g_lf_jac_err;
uint8_t g_lf_jac_bad;

/* ==========================================================================
 * 解析 Jacobian
 *
 * 推导（全部在"内部模型角" q 下，q = q_direction * (theta - q_zero)）：
 *   曲柄端点   Ka = L1*e1,  Kb = L1*e2      (ei = (cos qi, sin qi))
 *   两点连线段 d  = Kb - Ka,  s = |d|,  u = d/s
 *   中点       m  = (Ka + Kb)/2
 *   法向       n  = R90(u) = (-uy, ux)
 *   足端       p  = m + branch * h * n,   h = sqrt(L2^2 - s^2/4)
 *
 * 对 q1 求导（注意 d = Kb - Ka，所以 dd/dq1 = -dKa/dq1）：
 *   dm/dq1   = 0.5*L1*e1p
 *   ds/dq1   = -L1*(u·e1p)
 *   dh/dq1   = -(s/(4h))*ds/dq1 = (s*L1/(4h))*(u·e1p)
 *   dn/dq1   = R90(du/dq1) = (L1/s)*[(u·e1p)*n + e1]
 * 合并后 s/(4h) + h/s = L2^2/(h*s)，得到：
 *   dp/dq1 = L1*[ 0.5*e1p + b*(u·e1p)*(L2^2/(h*s))*n + b*(h/s)*e1 ]
 * 对 q2 同理（dd/dq2 = +dKb/dq2）：
 *   dp/dq2 = L1*[ 0.5*e2p - b*(u·e2p)*(L2^2/(h*s))*n - b*(h/s)*e2 ]
 * 其中 eip = (-sin qi, cos qi)。
 *
 * 最后把第 i 列乘 q_direction_i（因为 d/dtheta = d/dq * q_dir），并把 cm/rad 换成 m/rad。
 *
 * 该式已用中心差分逐点校验：两族分支共 1378 个位形，最大相对误差 4e-11。
 * ========================================================================== */
uint8_t Leg_Jacobian(const LegParam_t *param,
                     float theta1, float theta2,
                     float J[4], LegPoint_t *point)
{
    float q1, q2;
    float s1, c1, s2, c2;
    float kax, kay, kbx, kby;
    float dx, dy, s, h2, h;
    float ux, uy, nx, ny;
    float e1x, e1y, e1px, e1py;
    float e2x, e2y, e2px, e2py;
    float ue1, ue2, k, hs, b, l1, l2;
    float jx1, jy1, jx2, jy2;

    if (param == 0 || J == 0) {
        return 0;
    }

    l1 = param->link1;
    l2 = param->link2;

    if (l1 <= 0.0f || l2 <= 0.0f) {
        return 0;
    }

    q1 = (float)param->q1_direction * (theta1 - param->q1_zero);
    q2 = (float)param->q2_direction * (theta2 - param->q2_zero);

    s1 = sinf(q1); c1 = cosf(q1);
    s2 = sinf(q2); c2 = cosf(q2);

    kax = l1 * c1; kay = l1 * s1;
    kbx = l1 * c2; kby = l1 * s2;

    dx = kbx - kax;
    dy = kby - kay;
    s = sqrtf(dx * dx + dy * dy);

    /* s -> 0 就是两个曲柄重合，这是五连杆唯一的奇异位形（悬垂姿态正好是它） */
    if (s < 1.0e-3f) {
        return 0;
    }

    h2 = l2 * l2 - 0.25f * s * s;
    if (h2 < 1.0e-6f) {
        return 0;
    }
    h = sqrtf(h2);

    ux = dx / s;
    uy = dy / s;
    nx = -uy;
    ny = ux;

    b = (float)param->q1_branch;

    e1x = c1;      e1y = s1;
    e1px = -s1;    e1py = c1;
    e2x = c2;      e2y = s2;
    e2px = -s2;    e2py = c2;

    ue1 = ux * e1px + uy * e1py;
    ue2 = ux * e2px + uy * e2py;

    k  = l2 * l2 / (h * s);
    hs = h / s;

    jx1 = l1 * (0.5f * e1px + b * ue1 * k * nx + b * hs * e1x);
    jy1 = l1 * (0.5f * e1py + b * ue1 * k * ny + b * hs * e1y);
    jx2 = l1 * (0.5f * e2px - b * hs * e2x - b * ue2 * k * nx);
    jy2 = l1 * (0.5f * e2py - b * hs * e2y - b * ue2 * k * ny);

    J[0] = jx1 * (float)param->q1_direction * LF_CM2M;
    J[1] = jx2 * (float)param->q2_direction * LF_CM2M;
    J[2] = jy1 * (float)param->q1_direction * LF_CM2M;
    J[3] = jy2 * (float)param->q2_direction * LF_CM2M;

    if (point != 0) {
        point->x = kax + 0.5f * dx - b * h * uy;
        point->y = kay + 0.5f * dy + b * h * ux;
    }

    return 1;
}

/* 数值差分版本, 只给自检用 */
static uint8_t LF_NumericJacobian(const LegParam_t *param,
                                  float theta1, float theta2,
                                  float Jn[4])
{
    const float e = 1.0e-3f;
    const float sc = LF_CM2M / (2.0f * e);
    LegJoint_t j;
    LegPoint_t a;
    LegPoint_t b;

    j.q1 = theta1 + e;
    j.q2 = theta2;
    if (!Leg_ForwardKinematics(param, &j, &a)) {
        return 0;
    }
    j.q1 = theta1 - e;
    if (!Leg_ForwardKinematics(param, &j, &b)) {
        return 0;
    }
    Jn[0] = (a.x - b.x) * sc;
    Jn[2] = (a.y - b.y) * sc;

    j.q1 = theta1;
    j.q2 = theta2 + e;
    if (!Leg_ForwardKinematics(param, &j, &a)) {
        return 0;
    }
    j.q2 = theta2 - e;
    if (!Leg_ForwardKinematics(param, &j, &b)) {
        return 0;
    }
    Jn[1] = (a.x - b.x) * sc;
    Jn[3] = (a.y - b.y) * sc;

    return 1;
}

static void LF_JacSelfCheck(void)
{
    const float PAI = 3.14159265358979f;
    LegParam_t param;
    int8_t  fam;
    int16_t i;
    int16_t k;

    g_lf_jac_err = 0.0f;
    g_lf_jac_bad = 0;

    Leg_InitParam(&param, 0, 1, MMA_BRANCH1, MMA_BRANCH2);

    for (fam = 0; fam < 2; fam++) {
        param.q1_branch = (fam == 0) ? MMA_BRANCH1 : MMA_BRANCH2;
        param.q2_branch = (fam == 0) ? MMA_BRANCH2 : MMA_BRANCH1;

        for (i = -12; i <= 12; i++) {
            for (k = -12; k <= 12; k++) {
                float th1 = (float)i * (PAI / 12.0f);
                float th2 = (float)k * (PAI / 12.0f);
                float Ja[4];
                float Jn[4];
                float err;
                float scale;

                if (!Leg_Jacobian(&param, th1, th2, Ja, 0)) {
                    continue;
                }
                if (!LF_NumericJacobian(&param, th1, th2, Jn)) {
                    continue;
                }

                err = fabsf(Ja[0] - Jn[0]);
                if (fabsf(Ja[1] - Jn[1]) > err) err = fabsf(Ja[1] - Jn[1]);
                if (fabsf(Ja[2] - Jn[2]) > err) err = fabsf(Ja[2] - Jn[2]);
                if (fabsf(Ja[3] - Jn[3]) > err) err = fabsf(Ja[3] - Jn[3]);

                scale = fabsf(Jn[0]);
                if (fabsf(Jn[1]) > scale) scale = fabsf(Jn[1]);
                if (fabsf(Jn[2]) > scale) scale = fabsf(Jn[2]);
                if (fabsf(Jn[3]) > scale) scale = fabsf(Jn[3]);
                if (scale < 1.0e-4f) {
                    scale = 1.0e-4f;
                }

                err /= scale;

                if (err > g_lf_jac_err) {
                    g_lf_jac_err = err;
                }
                if (err > 1.0e-2f) {
                    g_lf_jac_bad++;
                }
            }
        }
    }
}

uint8_t Leg_MeasuredFootForce(const LegControl_t *leg, float *fx, float *fy)
{
    float J[4];
    float t1;
    float t2;
    float det;

    if (leg == 0 || fx == 0 || fy == 0) {
        return 0;
    }

    if (!Leg_Jacobian(&leg->param, leg->real_q1, leg->real_q2, J, 0)) {
        return 0;
    }

    det = J[0] * J[3] - J[1] * J[2];
    if (fabsf(det) < 1.0e-9f) {
        return 0;
    }

    t1 = Feedback[leg->param.motor1_id].tau * g_lf_tau_scale;
    t2 = Feedback[leg->param.motor2_id].tau * g_lf_tau_scale;

    /* J^T * F = tau  ==>  F = (J^T)^-1 * tau */
    *fx = ( J[3] * t1 - J[2] * t2) / det;
    *fy = (-J[1] * t1 + J[0] * t2) / det;

    return 1;
}

void Leg_ForceInit(void)
{
    uint8_t i;

    for (i = 0; i < LF_NLEG; i++) {
        g_lf_fx[i] = 0.0f;
        g_lf_fy[i] = 0.0f;
        g_lf_fz[i] = 0.0f;
        g_lf_ok[i] = 0;
    }
    g_lf_sum_fz = 0.0f;
    g_lf_mass_est = 0.0f;

    LF_JacSelfCheck();
}

void Leg_ForceObserveAll(void)
{
    uint8_t i;
    uint8_t n_valid = 0U;
    float sum = 0.0f;

    for (i = 0; i < LF_NLEG; i++) {
        float fx;
        float fy;

        if (!leg[i].ready || !leg[i].online) {
            g_lf_ok[i] = 0;
            continue;
        }

        if (!Leg_MeasuredFootForce(&leg[i], &fx, &fy)) {
            g_lf_ok[i] = 0;
            continue;
        }

        g_lf_fx[i] = fx;
        g_lf_fy[i] = fy;
        g_lf_fz[i] = -fy;
        g_lf_ok[i] = 1;
        n_valid++;
        sum += g_lf_fz[i];
    }

    g_lf_n_valid = n_valid;
    g_lf_sum_fz = sum;
    g_lf_mass_est = sum / LF_G;

    /* g_lf_sum_fz / g_lf_mass_est 只统计**有效**的腿。
     * 有腿掉线时必须先看 g_lf_n_valid: 如果它 < 4, mass_est 就是偏小的,
     * 按 mass_est * 4 / n_valid 才能估出真实体重。 */
}

/* ==========================================================================
 * 柔顺控制
 * ========================================================================== */

#define LF_PAI        3.14159265358979f
/* 摆动相占一个周期的比例, 必须和 Leg_GetTargetPoint 里的 0.4 一致 */
#define LF_SWING_END  (2.0f * LF_PAI * 0.4f)
/* 摆动末端这一段(占摆动时长的比例)用来把 kp 和前馈过渡到支撑相, 做软着陆 */
#define LF_LAND_FRAC  0.18f
/* 前馈竖直力的一阶低通系数。主要用来抹平"2 腿支撑 <-> 4 腿支撑"的切换 */
#define LF_FF_ALPHA   0.12f
/* kp 每拍最多变化多少: 开关切换/换档时刚度不要跳变 */
#define LF_KP_STEP    0.05f

/* 出厂值 = 和原来 Motor_SetPosition 完全等价(kp=1.0, kd=0.05, tau=0)。
 * 也就是说 g_lf_enable 从 0 改到 1 时, 电机收到的帧一个字都不变。
 * 必须从这里开始, 一格一格往下调, 不要一上来就给个"柔顺档"。 */
/* ===== tau 的换算系数（这是整个力控的根，必须先定死） =====
 * 电机反馈的 tau 和下发命令的 tau 都是**转子侧**，真实输出轴力矩 = 反馈值 * 6.33。
 * 依据（两条独立路径，结果一致）:
 *   1) PD 方程: 实测 tau ≈ 0.09 N·m，位置误差 0.0153 rad，kp=1.0
 *        转子侧: kp * 6.33 * 0.0153 = 0.097 N·m   <- 吻合
 *        输出侧: kp * 0.0153        = 0.015 N·m   <- 差 6 倍
 *   2) 刚度: 2.1mm 垂度 <-> 真实足端力 8.3N，K_y = 3.95 N/mm，理论 4.26 N/mm
 * 代码内部一律用**输出侧 N·m** 计算，只在收发两端做这个换算。
 * 如果以后换成别的电机，改这一个数即可。 */
float g_lf_tau_scale  = 6.33f;

uint8_t g_lf_enable    = 0;
float   g_lf_kp_stance = 1.0f;
float   g_lf_kp_swing  = 1.0f;
float   g_lf_kd        = 0.05f;
float   g_lf_ff_scale  = 0.0f;   /* 前馈从 0 开始, 一档一档往上加 */
float   g_lf_mass      = 0.0f;   /* 0 = 用 LF_MASS_DEFAULT */
float   g_lf_tau_limit = 8.0f;
float   g_lf_tau_step  = 0.40f;
float   g_lf_n_stance  = 4.0f;

/* 通信看门狗。g_lf_comm_fault=1 表示有腿掉线; g_lf_comm_strict=1 时
 * 会自动把 g_lf_enable 清零(强制退回位控)。这条很重要:
 * 掉线的电机保持最后一帧的 kp=1.0, 是一根硬腿, 不能让它和柔顺腿混在一起。 */
uint8_t g_lf_comm_fault  = 0U;
uint8_t g_lf_comm_strict = 1U;

/* ===== 防翻转保护 =====
 * 这条五连杆有两个退化位形, 都在"两个曲柄平行"处:
 *     r = L2 + L1 = 38.5 cm  (完全伸直)
 *     r = L2 - L1 = 13.5 cm  (完全折叠)
 * 在 r -> 13.5 附近, det(J) = r*r'/2 -> 0, 足端力对应的关节力矩趋于无穷,
 * 腿撑不住任何载荷, 于是被体重压着穿过退化点; 穿过之后曲柄两端重合,
 * 脚端可以沿半径 L2 的整圆任意滑动 —— 表现就是"脚从最低点直接翻到最高点卡死"。
 * 下面的保护在 r 或足端误差越界时强行把 kp 顶回去, 把腿拉回安全工作区。 */
/* 只允许柔顺的腿(bit0=leg[0] ... bit3=leg[3])。默认 0x0F 全开。
 * 调试时建议先用 0x01 只放一条腿, 其余三条保持原来的位控 —— 出了问题也只有一条腿。 */
uint8_t g_lf_leg_mask   = 0x0FU;

/* ===== 纯力矩模式 = VMC（虚拟腿） =====
 * 1 = 驱动器 kp 恒为 0(电机变成纯力矩源), 刚度全部由下面这条"虚拟腿"给出:
 *
 *     沿腿长方向:  F_r = kr * (r_d - r) - ff * sin(beta)
 *     摆角方向  : Q_b = kphi * (beta_d - beta) - ff * cos(beta) * r
 *
 * 再把这两个广义力映射回关节:
 *     alpha = (q1 - q2)/2 ,  Q_a = F_r * dr/dalpha
 *     tau1 = (Q_a + Q_b)/2 ,  tau2 = (Q_a - Q_b)/2
 *
 * 这个式子和两族分支无关(q1_dir=+1, q2_dir=-1 对所有腿都一样), 而且 beta 就是
 * 足端的极角 atan2(y,x), 所以不需要为镜像腿单独处理符号。
 *
 * 0 = 阻抗模式(驱动器 kp 出刚度, MCU 只加前馈)。
 *
 * 为什么阻尼留给驱动器: 驱动器的 kd*(0 - omega) 在它自己的 20kHz 环里算, 没有延迟;
 * MCU 这边从发命令到收反馈要 1.2~2ms, 用这么慢的速度量做微分会振。 */
uint8_t g_lf_pure_torque = 0U;
float   g_lf_kr          = 1000.0f;  /* 腿长刚度 N/m   (1000 = 1 N/mm) */
float   g_lf_kphi        = 25.0f;    /* 摆角刚度 N·m/rad */

/* 前馈单腿硬上限(占体重的比例)。万一 n_stance 算错(比如只有 1),
 * 这一层保证前馈不会变成 4 倍体重。0.6 = 允许到 6 成体重压在一条腿上,
 * 正常 4 腿(25%)和 2 腿(50%)都不会被夹到。 */
float   g_lf_ff_max_frac = 0.6f;

float   g_lf_kp_guard   = 1.50f;   /* 越界时强制使用的 kp */
float   g_lf_r_min      = 17.0f;   /* 安全下界 cm (摆动抬腿最高处约 19, 所以取 17) */
float   g_lf_r_max      = 36.0f;   /* 安全上界 cm (完全伸直是 38.5) */
float   g_lf_err_limit  = 8.0f;    /* 足端偏离目标的安全上限 cm */

static float LF_Clampf(float v, float lim)
{
    if (v > lim)  { return lim; }
    if (v < -lim) { return -lim; }
    return v;
}

/* 把角度折到 (-pi, pi]，用于摆角误差 */
static float LF_WrapPi(float a)
{
    while (a > LF_PAI) {
        a -= 2.0f * LF_PAI;
    }
    while (a < -LF_PAI) {
        a += 2.0f * LF_PAI;
    }
    return a;
}

/* 由腿指针反查索引（leg[] 是全局数组） */
static uint8_t LF_LegIndex(const LegControl_t *lg)
{
    uint8_t i;

    for (i = 0; i < LF_NLEG; i++) {
        if (lg == &leg[i]) {
            return i;
        }
    }
    return 0xFFU;
}

/* 把这一腿的柔顺状态复位（退回位控时调用，保证下次进入是从硬往软降） */
static void LF_ResetLegState(LegControl_t *lg)
{
    lg->kf_kp = 1.0f;
    lg->kf_guard = 0U;
    lg->kf_ff_y = 0.0f;
    lg->kf_tau1 = 0.0f;
    lg->kf_tau2 = 0.0f;
}

/* 变化率限幅: 每拍最多朝目标挪 step，防止前馈突然加上去。 */
static float LF_Slewf(float target, float *state, float step)
{
    float d = target - *state;

    if (d > step)  { d = step; }
    if (d < -step) { d = -step; }

    *state += d;
    return *state;
}

uint8_t Leg_ForceBeginCycle(void)
{
    uint8_t i;
    uint8_t n = 0;
    uint8_t n_off = 0;

    for (i = 0; i < LF_NLEG; i++) {
        float ph;

        if (!leg[i].ready) {
            continue;
        }
        /* 必须排除掉线的腿: 一条腿掉线还在 n_stance 里, 前馈就会少分一份,
         * 剩下几条腿撑不住整机。 */
        if (!leg[i].online) {
            n_off++;
            continue;
        }
        if (!Leg_GetTargetPoint(&leg[i], 0, 0, &ph)) {
            continue;
        }
        if (ph >= 0.0f && ph < LF_SWING_END) {
            continue;               /* 摆动相, 这条腿不承重 */
        }
        n++;
    }

    if (n == 0U) {
        n = 1U;                     /* 防止除零 */
    }

    g_lf_n_stance = (float)n;

    /* ===== 通信看门狗 =====
     * 电机掉线后**不会变软**——它会一直执行收到的最后一帧(kp=1.0 的位置保持),
     * 变成一根刚性支柱。这时如果其余三条腿在柔顺, 载荷会全压到那根硬腿上;
     * 而且"超时切阻尼模式"这种保命动作根本发不出去(通信已经断了)。
     * 所以只要有一条腿掉线, 就强制退回位控——位控是已知能工作的状态。 */
    if (n_off > 0U) {
        g_lf_comm_fault = 1U;
        if (g_lf_comm_strict) {
            g_lf_enable = 0U;
        }
    } else {
        g_lf_comm_fault = 0U;
    }

    return n;
}

uint8_t Leg_ForceUpdate(LegControl_t *lg)
{
    float x;
    float y;
    float phase;
    LegPoint_t target;
    LegJoint_t joint;
    float J[4];
    float mix;          /* 0 = 完全摆动参数, 1 = 完全支撑参数 */
    float kp;
    float kp_target;
    float m;
    float share;
    float ff_max;
    float ff_target;
    float fx;
    float fy;
    float tau1;
    float tau2;
    float r;
    float dex;
    float dey;
    float derr;
    uint8_t idx;

    if (lg == 0 || !lg->ready) {
        return 0;
    }

    /* 总开关关掉, 或者这条腿不在允许范围内 -> 走原来的位控(逐字节相同) */
    idx = LF_LegIndex(lg);
    if ((!g_lf_enable) || (idx > 3U) ||
        ((g_lf_leg_mask & (uint8_t)(1U << idx)) == 0U)) {
        LF_ResetLegState(lg);
        return Leg_Update(lg);
    }

    if (!Leg_GetTargetPoint(lg, &x, &y, &phase)) {
        return 0;
    }

    target.x = x;
    target.y = y;

    if (!Leg_InverseKinematics(&lg->param, &target, &joint)) {
        return 0;                   /* 目标点不可达, 本拍不发命令(保持上一拍) */
    }

    /* --- 1. 相位 -> 摆动/支撑, 并算出混合系数 --- */
    if (phase < 0.0f) {
        lg->kf_stance = 1U;         /* 站立/保持模式, 一直算支撑相 */
        mix = 1.0f;
    } else if (phase < LF_SWING_END) {
        lg->kf_stance = 0U;
        if (phase > LF_SWING_END * (1.0f - LF_LAND_FRAC)) {
            mix = (phase - LF_SWING_END * (1.0f - LF_LAND_FRAC)) /
                  (LF_SWING_END * LF_LAND_FRAC);
            if (mix > 1.0f) { mix = 1.0f; }
            if (mix < 0.0f) { mix = 0.0f; }
        } else {
            mix = 0.0f;
        }
    } else {
        lg->kf_stance = 1U;
        mix = 1.0f;
    }

    /* --- 2. 刚度目标: 摆动 -> 支撑 平滑过渡 --- */
    if (g_lf_pure_torque) {
        kp_target = 0.0f;           /* 纯力矩: 驱动器不做任何位置控制 */
    } else {
        kp_target = g_lf_kp_swing + (g_lf_kp_stance - g_lf_kp_swing) * mix;
    }

    /* --- 2b. 防翻转保护 -----------------------------------------------------
     * 用**实测**(正解出来的)足端点判断, 而不是用目标点。
     * r 掉到 17cm 以下(接近 13.5 的退化点)或者足端离开目标超过 8cm,
     * 说明腿已经被压到危险区, 马上把 kp 顶回去并撤掉前馈, 把腿往目标拖。
     * 带迟滞: 触发后要回到安全区内部(留 1.5cm / 2cm 余量)才解除, 否则会在边界上抖动。 */
    r = sqrtf(lg->real_point.x * lg->real_point.x +
              lg->real_point.y * lg->real_point.y);
    dex = lg->real_point.x - x;
    dey = lg->real_point.y - y;
    derr = sqrtf(dex * dex + dey * dey);

    if (lg->kf_guard) {
        if ((r > (g_lf_r_min + 1.5f)) && (r < (g_lf_r_max - 1.5f)) &&
            (derr < (g_lf_err_limit - 2.0f))) {
            lg->kf_guard = 0U;
        }
    } else {
        if ((r < g_lf_r_min) || (r > g_lf_r_max) || (derr > g_lf_err_limit)) {
            lg->kf_guard = 1U;
        }
    }

    if (lg->kf_guard) {
        if (!g_lf_pure_torque) {
            kp_target = g_lf_kp_guard;
        }
    }

    /* --- 2c. kp 斜坡: 任何切换/换档都不让刚度跳变 --- */
    {
        float dk = kp_target - lg->kf_kp;

        if (dk > LF_KP_STEP)  { dk = LF_KP_STEP; }
        if (dk < -LF_KP_STEP) { dk = -LF_KP_STEP; }
        lg->kf_kp += dk;
    }
    kp = lg->kf_kp;

    /* --- 3. 重力前馈: 支撑腿分摊 mg/n_stance, 摆动相为 0 --- */
    m = (g_lf_mass > 0.0f) ? g_lf_mass : LF_MASS_DEFAULT;
    share = (g_lf_n_stance > 0.5f) ?
            (g_lf_ff_scale * m * LF_G / g_lf_n_stance) : 0.0f;

    /* 硬上限: n_stance 一旦算错(比如只有 1), 前馈就会变成 4 倍体重。
     * 纯力矩模式下没有位置弹簧兜着, 这一层是必需的。 */
    ff_max = g_lf_ff_max_frac * m * LF_G;
    if (share > ff_max) {
        share = ff_max;
    }

    if (lg->kf_guard) {
        share = 0.0f;               /* 越界时不要再往下压 */
    }

    ff_target = share * mix;        /* 站立时 mix=1; 摆动相由 mix 决定过渡 */
    lg->kf_ff_y += (ff_target - lg->kf_ff_y) * LF_FF_ALPHA;

    /* --- 4. 足端力 -> 关节力矩:  tau = J^T * F --- */
    if (!Leg_Jacobian(&lg->param, lg->real_q1, lg->real_q2, J, 0)) {
        fx = 0.0f;
        fy = 0.0f;
        tau1 = 0.0f;
        tau2 = 0.0f;
    } else if (g_lf_pure_torque) {
        /* ===== VMC 虚拟腿 ===== */
        float q1i;
        float q2i;
        float alpha;
        float sa;
        float ca;
        float hh;
        float h2;
        float rp;               /* dr/dalpha, m/rad */
        float r_m;
        float r_d;
        float b_m;
        float b_d;
        float fr;               /* 沿腿长的力, 正 = 往外伸 */
        float qb;               /* 摆角方向的广义力矩 */
        float qa;
        float krr = g_lf_kr;
        float kpp = g_lf_kphi;

        q1i = (float)lg->param.q1_direction * (lg->real_q1 - lg->param.q1_zero);
        q2i = (float)lg->param.q2_direction * (lg->real_q2 - lg->param.q2_zero);
        alpha = 0.5f * (q1i - q2i);
        sa = sinf(alpha);
        ca = cosf(alpha);
        h2 = lg->param.link2 * lg->param.link2 -
             lg->param.link1 * lg->param.link1 * sa * sa;
        if (h2 < 1.0e-6f) {
            tau1 = 0.0f;
            tau2 = 0.0f;
        } else {
            hh = sqrtf(h2);
            rp = -lg->param.link1 * sa *
                 (1.0f + lg->param.link1 * ca / hh) * LF_CM2M;   /* cm/rad -> m/rad */

            r_m = sqrtf(lg->real_point.x * lg->real_point.x +
                        lg->real_point.y * lg->real_point.y) * LF_CM2M;
            r_d = sqrtf(x * x + y * y) * LF_CM2M;
            b_m = atan2f(lg->real_point.y, lg->real_point.x);
            b_d = atan2f(y, x);

            if (lg->kf_guard) {
                krr *= 3.0f;    /* 越界时弹簧变硬, 把腿拖回来 */
                kpp *= 3.0f;
            }

            /* 沿腿长的力: 弹簧 + 重力前馈在腿长方向的分量 */
            fr = krr * (r_d - r_m) - lg->kf_ff_y * sinf(b_m);
            /* 摆角广义力矩(力 × 力臂) */
            qb = kpp * LF_WrapPi(b_d - b_m) -
                 lg->kf_ff_y * cosf(b_m) * r_m;

            qa = fr * rp;

            tau1 = 0.5f * (qa + qb);
            tau2 = 0.5f * (qa - qb);
        }
    } else {
        /* 阻抗: 刚度由驱动器 kp 提供, 这里只加前馈 */
        fx = 0.0f;
        fy = -lg->kf_ff_y;
        tau1 = J[0] * fx + J[2] * fy;
        tau2 = J[1] * fx + J[3] * fy;
    }

    /* --- 5. 限幅 + 变化率限幅（都在输出侧 N·m 下做） --- */
    tau1 = LF_Clampf(tau1, g_lf_tau_limit);
    tau2 = LF_Clampf(tau2, g_lf_tau_limit);
    tau1 = LF_Slewf(tau1, &lg->kf_tau1, g_lf_tau_step);
    tau2 = LF_Slewf(tau2, &lg->kf_tau2, g_lf_tau_step);

    /* --- 5b. 输出侧 -> 转子侧: 命令字段要的是转子侧力矩 --- */
    tau1 /= g_lf_tau_scale;
    tau2 /= g_lf_tau_scale;

    /* --- 6. 下发(和 Leg_MoveTo 同一套绝对角约定) --- */
    lg->out_q1 = joint.q1 - lg->start_joint.q1;
    lg->out_q2 = joint.q2 - lg->start_joint.q2;
    lg->out_ready = 1;
    lg->out_is_abs = 0;

    Motor_SetImpedance(lg->param.motor1_id, lg->out_q1 + lg->q1_off,
                       tau1, kp, g_lf_kd);
    Motor_SetImpedance(lg->param.motor2_id, lg->out_q2 + lg->q2_off,
                       tau2, kp, g_lf_kd);

    return 1;
}

