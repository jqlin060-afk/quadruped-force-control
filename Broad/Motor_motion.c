#include "Motor_motion.h"
#include <math.h>

#define PAI_F 3.14159265358979323846f

LegControl_t leg[4];

/* ==========================================================================
 * 悬垂标定: 每条腿实际的"劈叉角" alpha (rad)。完全伸直 = 0。
 *
 * 标定时的假设是"四条腿都完全伸直"。实际上腿可能被线束、机械限位挡住,
 * 只伸到 r < 38.5。这时如果按 0 来算, 那条腿的电机零位就偏了, 而且:
 *
 *   r 在 alpha=0 附近对角度极不敏感  ->  1 cm 的长度差对应约 19 度!
 *
 * 查表(尺子量到"脚比最长那条低多少"):
 *     低 0.2 cm -> 0.147      低 1.0 cm -> 0.331      低 3.0 cm -> 0.584
 *     低 0.5 cm -> 0.233      低 1.5 cm -> 0.408      低 4.0 cm -> 0.680
 *     低 0.8 cm -> 0.300      低 2.0 cm -> 0.473      低 6.0 cm -> 0.850
 *
 * 量法: 狗悬空、腿自然垂到底, 拿一把直尺横在四个脚底, 量每只脚到同一基准
 * (比如机身下表面)的高度。取最长的那条为基准(r=38.5), 其余的查表填这里。
 * ========================================================================== */
static const float LF_HANG_ALPHA[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

/* 四条腿的模型参数。起始姿态给"悬垂满伸", 因为上电时腿是自然垂直朝下的 */
void Led_ControlUnitInit()
{
  Leg_ControlInit(&leg[0], FL_M1, FL_M2, 1,  -1, 0.0f, MMA_HANG_Y);  /* 左前 */
  Leg_ControlInit(&leg[1], FR_M1, FR_M2, -1,  1, 0.0f, MMA_HANG_Y);  /* 右前 */
  Leg_ControlInit(&leg[2], RL_M1, RL_M2, -1,  1, 0.0f, MMA_HANG_Y);  /* 左后 */
  Leg_ControlInit(&leg[3], RR_M1, RR_M2, 1,  -1, 0.0f, MMA_HANG_Y);  /* 右后 */

  /* 逐腿的悬垂劈叉角, 见上面 LF_HANG_ALPHA 的说明 */
  leg[0].param.hang_alpha = LF_HANG_ALPHA[0];
  leg[1].param.hang_alpha = LF_HANG_ALPHA[1];
  leg[2].param.hang_alpha = LF_HANG_ALPHA[2];
  leg[3].param.hang_alpha = LF_HANG_ALPHA[3];
}

static float ClampUnit(float v)
{
    if (v > 1.0f) {
        return 1.0f;
    }
    if (v < -1.0f) {
        return -1.0f;
    }
    return v;
}

uint8_t Leg_ForwardKinematics(const LegParam_t *param,
                              const LegJoint_t *joint,
                              LegPoint_t *point)
{
    float q1;
    float q2;
    float kax;
    float kay;
    float kbx;
    float kby;
    float dx;
    float dy;
    float d;
    float half;
    float h2;
    float h;
    float ux;
    float uy;
    float mx;
    float my;
    float branch;

    if (param == 0 || joint == 0 || point == 0) {
        return 0;
    }

    if (param->link1 <= 0.0f || param->link2 <= 0.0f) {
        return 0;
    }

    q1 = (float)param->q1_direction * (joint->q1 - param->q1_zero);
    q2 = (float)param->q2_direction * (joint->q2 - param->q2_zero);

    kax = param->link1 * cosf(q1);
    kay = param->link1 * sinf(q1);
    kbx = param->link1 * cosf(q2);
    kby = param->link1 * sinf(q2);

    dx = kbx - kax;
    dy = kby - kay;
    d = sqrtf(dx * dx + dy * dy);

    if (d < MMA_EPS) {
        return 0;
    }

    if (d > 2.0f * param->link2 + MMA_EPS) {
        return 0;
    }

    half = 0.5f * d;
    h2 = param->link2 * param->link2 - half * half;

    if (h2 < 0.0f) {
        h2 = 0.0f;
    }

    h = sqrtf(h2);
    ux = dx / d;
    uy = dy / d;
    mx = kax + half * ux;
    my = kay + half * uy;

    branch = (float)param->q1_branch;

    point->x = mx - branch * h * uy;
    point->y = my + branch * h * ux;

    return 1;
}

uint8_t Leg_InverseKinematics(const LegParam_t *param,
                              const LegPoint_t *point,
                              LegJoint_t *joint)
{
    float r;
    float alpha;
    float base;
    float cosine;

    if (param == 0 || point == 0 || joint == 0) {
        return 0;
    }

    if (param->link1 <= 0.0f || param->link2 <= 0.0f) {
        return 0;
    }

    r = sqrtf(point->x * point->x + point->y * point->y);

    if (r < MMA_EPS) {
        return 0;
    }

    if (r > param->link1 + param->link2 + MMA_EPS ||
        r < fabsf(param->link1 - param->link2) - MMA_EPS) {
        return 0;
    }

    cosine = (param->link1 * param->link1 + r * r -
              param->link2 * param->link2) /
             (2.0f * param->link1 * r);
    alpha = acosf(ClampUnit(cosine));
    base = atan2f(point->y, point->x);

    joint->q1 = param->q1_zero +
                (float)param->q1_direction *
                (base + (float)param->q1_branch * alpha);

    joint->q2 = param->q2_zero +
                (float)param->q2_direction *
                (base + (float)param->q2_branch * alpha);

    return 1;
}

void Leg_InitParam(LegParam_t *param,
                   uint8_t motor1_id,
                   uint8_t motor2_id,
                   int8_t q1_branch,
                   int8_t q2_branch)
{
    if (param == 0) {
        return;
    }

    param->link1 = MMA_L1;
    param->link2 = MMA_L2;

    param->q1_zero = 0.0f;
    param->q2_zero = 0.0f;
    param->q1_direction = 1;
    param->q2_direction = -1;

    param->q1_branch = q1_branch;
    param->q2_branch = q2_branch;

    param->hang_alpha = MMA_HANG_ALPHA;    /* 默认 0, 由 Led_ControlUnitInit 逐腿覆盖 */

    param->motor1_id = motor1_id;
    param->motor2_id = motor2_id;
}

static void Leg_SendBoth(uint8_t motor1_id,
                         float theta1,
                         uint8_t motor2_id,
                         float theta2)
{
    Motor_SetPosition(motor1_id, theta1);
    Motor_SendCmdDMA();

    Motor_SetPosition(motor2_id, theta2);
    Motor_SendCmdDMA();
}

uint8_t Leg_ControlInit(LegControl_t *leg,
                        uint8_t motor1_id,
                        uint8_t motor2_id,
                        int8_t q1_branch,
                        int8_t q2_branch,
                        float start_x,
                        float start_y)
{
    LegPoint_t start_point;

    if (leg == 0) {
        return 0;
    }

    leg->ready = 0;
    Leg_InitParam(&leg->param, motor1_id, motor2_id, q1_branch, q2_branch);

    start_point.x = start_x;
    start_point.y = start_y;

    if (!Leg_InverseKinematics(&leg->param,
                               &start_point,
                               &leg->start_joint)) {
        return 0;
    }

    Motor_ResetPositionControl(motor1_id);
    Motor_ResetPositionControl(motor2_id);

    leg->q1_off = 0.0f;
    leg->q2_off = 0.0f;

    leg->out_ready = 0;
    leg->out_is_abs = 0;
    leg->out_q1 = 0.0f;
    leg->out_q2 = 0.0f;

    leg->real_q1 = 0.0f;
    leg->real_q2 = 0.0f;
    leg->real_point.x = start_x;
    leg->real_point.y = start_y;

    leg->err_q1 = 0.0f;
    leg->err_q2 = 0.0f;
    leg->online = 0;

    leg->kf_tau1 = 0.0f;
    leg->kf_tau2 = 0.0f;
    leg->kf_ff_y = 0.0f;
    leg->kf_stance = 0;
    leg->kf_contact = 0;
    leg->kf_guard = 0;
    leg->kf_primed = 0;

    leg->ready = 1;

    return 1;
}

/* 上电悬垂标定。
 * 悬垂时足端在电机轴正下方, 两个曲柄关于竖直线对称: q1 = -90 度 + alpha, q2 = -90 度 - alpha。
 * 由 q = q_direction * (theta - q_zero) 反推零位:
 *     q_zero = theta - q / q_direction          (q_direction = +-1)
 * 精度只取决于 MMA_HANG_ALPHA 这个假设(见头文件里的说明), 跟这条腿前面的历史无关。 */
uint8_t Leg_HangCalib(LegControl_t *leg)
{
    float q1_hang;
    float q2_hang;

    if (leg == 0 || !leg->ready) {
        return 0;
    }

    Leg_ReadFeedback(leg);          /* 拿最新的 theta 和 online */

    if (!leg->online) {
        return 0;                   /* 电机还没上报, 角度不可信 */
    }

    /* 悬垂时两条大腿关于竖直线对称, 但对称方向由 branch 决定:
     *   branch (+1,-1): 左大腿 = -90 + alpha,  右大腿 = -90 - alpha
     *   branch (-1,+1): 左大腿 = -90 - alpha,  右大腿 = -90 + alpha
     * 必须按 branch 分开。否则 (-1,+1) 那两条腿(右前/左后)会多收 2*alpha,
     * 收起来比另外两条高, 站立就歪七扭八。 */
    q1_hang = -0.5f * PAI_F + (float)leg->param.q1_branch * leg->param.hang_alpha;
    q2_hang = -0.5f * PAI_F + (float)leg->param.q2_branch * leg->param.hang_alpha;

    leg->param.q1_zero = leg->real_q1 -
                         q1_hang / (float)leg->param.q1_direction;
    leg->param.q2_zero = leg->real_q2 -
                         q2_hang / (float)leg->param.q2_direction;

    /* 当前姿态就是运动基准。注意: 电机的 theta 命令是"绝对目标位置",
     * 不是"相对增量"。所以把 q_off 设成当前读到的 theta,
     * 这样 下发值 = (IK - start_joint) + q_off = IK, 就是绝对角。
     * 如果 q_off 留 0, 那"保持当前姿态"这一帧会发 0, 等于命令电机转到绝对 0,
     * 每个电机都会转掉自己上电时的那个 theta -> 四条腿各偏一个不同的量 -> 歪七扭八。 */
    leg->start_joint.q1 = leg->param.q1_zero +
                          q1_hang / (float)leg->param.q1_direction;
    leg->start_joint.q2 = leg->param.q2_zero +
                          q2_hang / (float)leg->param.q2_direction;

    leg->q1_off = leg->start_joint.q1;
    leg->q2_off = leg->start_joint.q2;

    /* 平滑器也要对齐到当前真实角度, 否则第一段会从 0 斜坡上来 */
    Motor_SyncPosition(leg->param.motor1_id, leg->real_q1);
    Motor_SyncPosition(leg->param.motor2_id, leg->real_q2);

    leg->out_ready = 1;
    leg->out_is_abs = 0;
    leg->out_q1 = 0.0f;
    leg->out_q2 = 0.0f;

    return 1;
}

void Leg_ReadFeedback(LegControl_t *leg)
{
    LegJoint_t joint;
    LegJoint_t abs_joint;

    if (leg == 0 || !leg->ready) {
        return;
    }

    leg->real_q1 = Feedback[leg->param.motor1_id].theta;
    leg->real_q2 = Feedback[leg->param.motor2_id].theta;

    leg->online = (uint8_t)(Motor_IsOnline(leg->param.motor1_id) &&
                            Motor_IsOnline(leg->param.motor2_id));

    /* 电机 theta 就是绝对关节角(相对上电零点)。
     * 对齐时下发的是绝对角, 平时下发的是相对基准的角。 */
    if (leg->out_is_abs) {
        leg->err_q1 = leg->out_q1 - leg->real_q1;
        leg->err_q2 = leg->out_q2 - leg->real_q2;
    } else {
        leg->err_q1 = leg->out_q1 - (leg->real_q1 - leg->start_joint.q1);
        leg->err_q2 = leg->out_q2 - (leg->real_q2 - leg->start_joint.q2);
    }

    abs_joint.q1 = leg->real_q1;
    abs_joint.q2 = leg->real_q2;

    joint = abs_joint;
    Leg_ForwardKinematics(&leg->param, &joint, &leg->real_point);
}

/* 读四条腿的反馈, 返回在线的腿数(0..4) */
uint8_t Leg_ReadAll(void)
{
    uint8_t i;
    uint8_t online_num = 0U;

    for (i = 0U; i < 4U; i++) {
        Leg_ReadFeedback(&leg[i]);
        if (leg[i].online) {
            online_num++;
        }
    }

    return online_num;
}

uint8_t Leg_MoveTo(LegControl_t *leg, float x, float y)
{
    LegPoint_t target;
    LegJoint_t joint;

    if (leg == 0 || !leg->ready) {
        return 0;
    }

    target.x = x;
    target.y = y;

    if (!Leg_InverseKinematics(&leg->param, &target, &joint)) {
        return 0;
    }

    leg->out_q1 = joint.q1 - leg->start_joint.q1;
    leg->out_q2 = joint.q2 - leg->start_joint.q2;
    leg->out_ready = 1;
    leg->out_is_abs = 0;

    Leg_SendBoth(leg->param.motor1_id, leg->out_q1 + leg->q1_off,
                 leg->param.motor2_id, leg->out_q2 + leg->q2_off);

    return 1;
}

uint8_t Leg_GetTargetPoint(LegControl_t *leg, float *x, float *y, float *phase)
{
    float t;
    float ph;
    float px;
    float py;

    if (leg == 0 || !leg->ready) {
        return 0;
    }

    if (leg->mode != LEG_MODE_SWING) {
        if (x != 0) { *x = leg->hold_x; }
        if (y != 0) { *y = leg->hold_y; }
        if (phase != 0) { *phase = -1.0f; }
        return 1;
    }

    t = (float)(HAL_GetTick() - leg->swing_t0) * 0.001f;
    ph = 2.0f * PAI_F * leg->swing_freq * t + leg->swing_phase;

    while (ph >= 2.0f * PAI_F) {
        ph -= 2.0f * PAI_F;
    }
    while (ph < 0.0f) {
        ph += 2.0f * PAI_F;
    }

    if (ph < (2.0f * PAI_F * 0.4f)) {
        float u = ph / (2.0f * PAI_F * 0.4f);
        float w = 2.0f * PAI_F * u;

        px = leg->swing_cx - leg->swing_ax +
             leg->swing_ax * (w - sinf(w)) / PAI_F;
        py = leg->swing_cy + leg->swing_ay * 0.5f * (1.0f - cosf(w));
    } else {
        float u = (ph - 2.0f * PAI_F * 0.4f) / (2.0f * PAI_F * 0.6f);

        px = leg->swing_cx + leg->swing_ax - 2.0f * leg->swing_ax * u;
        py = leg->swing_cy;
    }

    if (x != 0) { *x = px; }
    if (y != 0) { *y = py; }
    if (phase != 0) { *phase = ph; }

    return 1;
}

uint8_t Leg_Update(LegControl_t *leg)
{
    float x;
    float y;

    if (!Leg_GetTargetPoint(leg, &x, &y, 0)) {
        return 0;
    }

    return Leg_MoveTo(leg, x, y);
}
