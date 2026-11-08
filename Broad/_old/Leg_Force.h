#ifndef LEG_FORCE_H
#define LEG_FORCE_H

#include "Motor_motion.h"

/* ==========================================================================
 * 力控第一步: 纯观测层
 *
 * 这一层只"看"，不改任何控制量:
 *   - 五连杆解析 Jacobian（已用数值差分逐点校验）
 *   - 由关节力矩反馈 tau 反解足端力 F = (J^T)^-1 * tau
 *   - 四腿竖直力之和 / 机体重量估计
 *   - 反馈帧里原本被丢掉的 MERROR / FORCE 12bit 也一起解出来了（见 Motor.c）
 *
 * 悬空测试时能验证的东西:
 *   1) Jacobian 的分支符号对不对（加前馈时四条腿是同向伸还是两伸两缩）
 *   2) tau 的零点、噪声、符号
 *   3) MERROR / FORCE 字段到底动不动
 * 悬空测试**不能**验证: 机体重量的绝对值、足端刚度、抗推、上坡。
 * ========================================================================== */

#define LF_G             9.81f
/* 机体重量默认值(kg)。20~30 斤取中值 12.5；站稳之后的实测值会覆盖它。 */
#define LF_MASS_DEFAULT  12.5f

#define LF_NLEG          4u
#define LF_CM2M          0.01f     /* 运动学用的是 cm，力和力矩换算要 m */

/* ===== 观测量（调试器 Watch 里直接加这几个） ===== */
extern float   g_lf_fx[LF_NLEG];   /* 腿对地面的力: x 前后方向, N */
extern float   g_lf_fy[LF_NLEG];   /* 腿对地面的力: y 竖直方向, 承重时 < 0 (模型 y 向下为负) */
extern float   g_lf_fz[LF_NLEG];   /* = -g_lf_fy, 承重时 > 0, 看这个最直观 */
extern float   g_lf_sum_fz;        /* 四腿竖直力之和; 站在地上时应该 ≈ m*g ≈ 123 N */
extern float   g_lf_mass_est;      /* = g_lf_sum_fz / LF_G; 站在地上时就是体重的实测值 */
extern uint8_t g_lf_n_valid;       /* 上面两个数只统计有效的腿; < 4 时 mass_est 偏小,
                                    * 真实体重 ≈ g_lf_mass_est * 4 / g_lf_n_valid */
extern uint8_t g_lf_ok[LF_NLEG];   /* 0 = 该腿 J 奇异或离线(力不可信), 1 = 有效 */
extern float   g_lf_jac_err;       /* 解析 J vs 数值 J 的最大相对误差, 正常应 < 1e-2 */
extern uint8_t g_lf_jac_bad;       /* 自检失败的位形个数 */

/* 五连杆解析 Jacobian。
 * 输入 theta1/theta2 是**电机角**(和你下发的量同一套, 即 LegControl_t::real_q1/real_q2)。
 * 输出 J 行主序: J[0]=dx/dth1, J[1]=dx/dth2, J[2]=dy/dth1, J[3]=dy/dth2, 单位 m/rad。
 * point 可为 NULL; 非 NULL 时顺便把正解点(单位 cm)写出来。
 * 返回 0 = 位形无效(奇异), 此时 J 不可用。 */
uint8_t Leg_Jacobian(const LegParam_t *param,
                     float theta1, float theta2,
                     float J[4], LegPoint_t *point);

/* 由反馈的关节力矩反解足端力: J^T * F = tau  ==>  F = (J^T)^-1 * tau
 * (fx, fy) 是**腿对地面**的力, 单位 N。承重时 fy < 0。
 * 返回 0 = 该腿 J 奇异或电机离线, 输出未改动。 */
uint8_t Leg_MeasuredFootForce(const LegControl_t *leg, float *fx, float *fy);

/* 上电调用一次: 跑 Jacobian 自检（结果看 g_lf_jac_err / g_lf_jac_bad）。 */
void Leg_ForceInit(void);

/* 主循环里每拍调用一次: 刷新上面那组观测量。**不做任何控制**。 */
void Leg_ForceObserveAll(void);

/* ==========================================================================
 * 柔顺控制（第一步: 驱动器内阻抗 + 重力前馈）
 *
 * 控制律完全在驱动器里闭合（电机手册 §7.7 力位混合模式）:
 *     tau_out = tau_ff + kp*(theta_des - theta) + kd*(0 - omega)
 * 我们只负责每拍给出 theta_des(=IK(p_d))、kp/kd(按相位切换)、tau_ff(=J^T*F_ff)。
 * 好处: 稳定裕度还是由驱动器 20kHz 的环保证，MCU 1~2ms 的周期只是"缓慢改弹簧和偏置"，
 *       几乎不可能振起来；而 K 一降下来，手感/抗冲击能力立刻变成"柔顺"。
 *
 * 调参顺序（**必须按这个顺序，一步一档**）:
 *
 *   0) 先标定刚度尺度（不做这一步，任何 kp 都是瞎猜）:
 *      保持 g_lf_enable = 0（原位控），狗站地上不动，读 leg[i].real_point.y。
 *        ≈ -26.x           -> 足端竖直刚度约 4.2 N/mm（K_joint = 40*kp）
 *        ≈ -30 / -31       -> 足端竖直刚度约 0.67 N/mm（K_joint = 6.33*kp，差 6.33 倍）
 *      下沉量的绝对值直接告诉你"现在到底有多硬"，也顺便定掉了 tau 的换算系数。
 *
 *   1) 把 g_lf_enable 改成 1。**此时行为应该和 0 时一模一样**（kp 出厂就是 1.0，
 *      前馈是 0），腿不动、不塌、不抖。这一步只是在验证管路通了。
 *
 *   2) 一格一格降 g_lf_kp_stance（和 g_lf_kp_swing 一起）: 1.0 -> 0.9 -> 0.8 -> 0.7 ...
 *      每降一格, 看 real_point.y 往下掉了多少、kf_guard 有没有变 1。
 *      kf_guard 一变 1 就说明已经到安全边界了，退回上一档。
 *      想要"站着不塌但更软"，就要配合第 3 步的前馈一起加。
 *
 *   3) g_lf_ff_scale 从 0 开始: 0.1 -> 0.2 -> 0.3 ... 一档一档加。
 *      每档看 real_point.y 回到 STAND_Y 了没有（回来=前馈够了）。
 *
 *   4) 走路: 先 0.8Hz 原步态, 观察腿有没有拖地/抖。
 * ========================================================================== */

/* 运行总开关。0 = 完全退回原来的位控(Leg_Update)，行为和现在一模一样。
 * 调试器里直接改这个变量就能切换, 不用重新编译。 */
extern uint8_t g_lf_enable;

/* 只允许柔顺的腿: bit0=leg[0] ... bit3=leg[3]。
 * **第一次测试请用 0x01**(只放左前一条腿), 其余三条保持原来的位控。
 * 这样即使参数不对, 也只有一条腿会动, 整机不会塌。 */
extern uint8_t g_lf_leg_mask;

/* 支撑相 / 摆动相的电机 kp。**出厂值 1.0 = 完全不改现有行为**，
 * 这个数就是"软硬旋钮"：越小越柔顺。
 * 一格一格往下调（0.9 -> 0.8 -> ...），每降一格看 leg[i].real_point.y 有没有往下掉。 */
extern float g_lf_kp_stance;
extern float g_lf_kp_swing;
extern float g_lf_kd;

/* 重力前馈缩放: 0 = 不发力(纯柔顺, 站着会下沉), 1.0 = 按体重理论值发力。
 * **必须从 0 开始往上加。** */
extern float g_lf_ff_scale;

/* 机体质量(kg)。0 = 用 LF_MASS_DEFAULT。
 * 站到地上之后 g_lf_mass_est 就是实测体重, 抄过来即可。 */
extern float g_lf_mass;

/* tau 的换算系数: 真实输出轴力矩 = 电机反馈值 × 这个数。
 * 已用 PD 方程 + 足端刚度两条独立路径确认为 6.33（转子侧 / 减速比）。
 * 代码内部一律用输出侧 N·m 计算，收发两端各自换算。 */
extern float g_lf_tau_scale;

/* ===== VMC 纯力矩模式（虚拟腿） =====
 * 1 = 驱动器 kp 恒为 0(电机变成纯力矩源), 刚度由"虚拟腿"给出:
 *       沿腿长: F_r = kr*(r_d - r) + 重力分量
 *       摆角  : Q_b = kphi*(beta_d - beta) + 重力分量
 *       tau1 = (F_r*dr/da + Q_b)/2 ,  tau2 = (F_r*dr/da - Q_b)/2
 * 0 = 阻抗模式(驱动器 kp 出刚度, MCU 只加前馈)。
 *
 * 阻尼仍然交给驱动器的 kd —— 它在 20kHz 环里算, 没有通信延迟;
 * MCU 这边从发命令到收反馈要 1.2~2ms, 用这么慢的速度量做微分会振。 */
extern uint8_t g_lf_pure_torque;
extern float   g_lf_kr;            /* 腿长刚度 N/m,  1000 = 1 N/mm */
extern float   g_lf_kphi;          /* 摆角刚度 N·m/rad */

/* 前馈单腿硬上限(占体重比例)。防住 n_stance 算错时前馈变成 4 倍体重。 */
extern float   g_lf_ff_max_frac;

/* 前馈力矩的安全限幅 */
extern float g_lf_tau_limit;   /* 单电机 Nm */
extern float g_lf_tau_step;    /* 每拍最大变化量 Nm */

/* ===== 防翻转保护 =====
 * 五连杆在 r = 13.5cm(完全折叠) 和 r = 38.5cm(完全伸直) 两处退化，
 * r -> 13.5 附近 det(J) -> 0，腿撑不住载荷会被压过退化点，
 * 之后脚端能沿半径 26cm 的整圆乱滑 —— 就是"脚翻到最高点卡死"。
 * 下面四个参数定义安全工作区，越界就强行把 kp 顶到 g_lf_kp_guard。 */
extern float g_lf_kp_guard;    /* 越界时强制使用的 kp，默认 1.5 */
extern float g_lf_r_min;       /* 安全下界 cm，默认 17 (摆动抬腿最高约 19) */
extern float g_lf_r_max;       /* 安全上界 cm，默认 36 (完全伸直 38.5) */
extern float g_lf_err_limit;   /* 足端偏离目标的安全上限 cm，默认 7 */

/* 当前在支撑相的腿数(用于分摊体重), 调试用 */
extern float g_lf_n_stance;

/* ===== 通信看门狗 =====
 * 电机掉线后不会变软: 它一直执行最后一帧(kp=1.0 的位置保持), 变成一根刚性支柱,
 * 而其余腿在柔顺 -> 载荷全压到那根硬腿上; 同时"超时切阻尼模式"也发不出去。
 * 所以只要有一条腿掉线, 就强制退回位控。 */
extern uint8_t g_lf_comm_fault;    /* 1 = 有腿掉线(只读) */
extern uint8_t g_lf_comm_strict;   /* 1 = 掉线时自动把 g_lf_enable 清零(默认) */

/* 每拍在 4 条腿更新之前调用一次: 统计当前有几条腿在支撑相。 */
uint8_t Leg_ForceBeginCycle(void);

/* 柔顺版 Leg_Update。g_lf_enable=0 时内部直接转调 Leg_Update。 */
uint8_t Leg_ForceUpdate(LegControl_t *leg);

#endif
