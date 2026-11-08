#ifndef _MOTOR_MOTION_H_
#define _MOTOR_MOTION_H_

#include "Motor_advanced.h"
#include "Motor.h"

#define MMA_L1      12.5f

#define MMA_L2      26.0f

#define MMA_R_MAX   38.0f
#define MMA_R_MIN   15.2f

#define MMA_BRANCH1  1
#define MMA_BRANCH2 (-1)

#define MMA_EPS     0.0001f

//电机ID
#define FL_M1    2
#define FL_M2    3
#define FR_M1    4
#define FR_M2    5
#define RL_M1    9
#define RL_M2    8
#define RR_M1    7
#define RR_M2    6

//站立高度
#define STAND_Y        (-26.0f)

/* 悬垂(上电)姿态: 狗悬空、腿自然垂直朝下, 足端正好在电机轴正下方。
 * 完全伸直时两个曲柄都竖直向下(q1 = q2 = -90 度), 足端就在这个 y 上。 */
#define MMA_HANG_Y     (-(MMA_L1 + MMA_L2))

/* 悬垂时两个曲柄的"劈叉角"alpha(rad): q1 = -90 度 + alpha, q2 = -90 度 - alpha。
 * alpha = 0 表示腿完全伸直(r = L1 + L2 = 38.5)。
 * 注意: 满伸附近曲柄角对伸长量极敏感(那里 dq/dr 趋于无穷, 实测 r 每差 1mm 角度约差 0.47rad),
 * 所以如果腿被机械限位或者线束挡着没伸到底, 要把按实测 r 反算出来的 alpha 填在这里:
 *     alpha = acos((L1^2 + r^2 - L2^2) / (2 * L1 * r))  */
#define MMA_HANG_ALPHA  0.0f


typedef struct
{
    float x;
    float y;
} LegPoint_t;

typedef struct
{

    float q1;
    float q2;
} LegJoint_t;

typedef struct
{

    float link1;

    float link2;

    float q1_zero;
    float q2_zero;
    int8_t q1_direction;
    int8_t q2_direction;

    int8_t q1_branch;
    int8_t q2_branch;

    uint8_t motor1_id;
    uint8_t motor2_id;
} LegParam_t;

typedef struct
{
    LegParam_t param;
    LegJoint_t start_joint;

    uint8_t mode;

    float swing_cx;
    float swing_cy;
    float swing_ax;
    float swing_ay;
    float swing_freq;
    float swing_phase;

    uint32_t swing_t0;

    float hold_x;
    float hold_y;

    /* 直接角度偏移(rad): 站立和运动时都加在下发的两个 theta 上。
     * 绕过运动学, 直接指定物理姿态。
     * 两个电机同时加约 ±2.954 (169.24 度) 可翻转关节方向。 */
    float q1_off;
    float q2_off;

    uint8_t out_ready;
    uint8_t out_is_abs;
    float out_q1;
    float out_q2;

    float real_q1;
    float real_q2;
    LegPoint_t real_point;

    float err_q1;
    float err_q2;

    uint8_t online;

    uint8_t ready;

    /* ===== 力控/柔顺状态（由 Leg_Force.c 维护，方便在调试器里逐腿观察） ===== */
    float   kf_tau1;        /* 上一拍的前馈力矩（输出侧 N·m，下发时再除以 g_lf_tau_scale） */
    float   kf_tau2;
    float   kf_ff_y;        /* 低通后的前馈竖直力, 正 = 腿在承重 (N) */
    float   kf_kp;          /* 当前实际使用的电机 kp（斜坡后的值） */
    uint8_t kf_stance;      /* 1 = 支撑相, 0 = 摆动相 */
    uint8_t kf_guard;       /* 1 = 触发了防翻转保护（被压实/撞限位） */
} LegControl_t;

/* 四条腿的控制对象, 定义在 Motor_motion.c */
extern LegControl_t leg[4];

#define LEG_MODE_HOLD    0U
#define LEG_MODE_SWING   1U

uint8_t Leg_AlignPush(LegControl_t *leg, float x, float y);

void Leg_AlignRelease(LegControl_t *leg, float x, float y);

void Leg_SetBase(LegControl_t *leg, float x, float y);

typedef struct
{
    uint8_t tested;

    uint8_t passed;

    float   max_err;

    float q1_in[8];

    float q2_in[8];
    float q1_out[8];

    float q2_out[8];
    float foot_x[8];

    float foot_y[8];
    float err[8];

} LegCheck_t;

void Leg_InitParam(LegParam_t *param,
                   uint8_t motor1_id,
                   uint8_t motor2_id,
                   int8_t q1_branch,
                   int8_t q2_branch);

void Led_ControlUnitInit(void);

uint8_t Leg_ForwardKinematics(const LegParam_t *param,
                              const LegJoint_t *joint,
                              LegPoint_t *point);

uint8_t Leg_InverseKinematics(const LegParam_t *param,
                              const LegPoint_t *point,
                              LegJoint_t *joint);

uint8_t Leg_ControlInit(LegControl_t *leg,
                        uint8_t motor1_id,
                        uint8_t motor2_id,
                        int8_t q1_branch,
                        int8_t q2_branch,
                        float start_x,
                        float start_y);

uint8_t Leg_MoveTo(LegControl_t *leg, float x, float y);

uint8_t Leg_Update(LegControl_t *leg);

/* 只算足端目标点, 不下发任何命令。
 * 输出 x/y 单位 cm(与 LegPoint_t 一致); phase 是 [-0, 2pi) 的步态相位,
 * 返回 -1 表示该腿不在摆动模式(站立/保持), 此时 x/y 取 hold_x/hold_y。
 * 任何输出指针都可以传 NULL。 */
uint8_t Leg_GetTargetPoint(LegControl_t *leg, float *x, float *y, float *phase);

void Leg_ReadFeedback(LegControl_t *leg);

uint8_t Leg_CheckKinematics(LegCheck_t *report);

/* 上电悬垂标定: 把"当前物理姿态"声明成模型里的悬垂姿态(腿垂直朝下、满伸),
 * 由此定两个电机的零位(theta = 0 时曲柄指向哪), 并把相对运动的基准 start_joint 挪到当前姿态。
 * 必须在腿真的悬垂、且已经读到电机反馈之后调用。
 * 返回 1 = 标定完成; 0 = 腿未就绪或电机未上线(此时不改动任何参数)。 */
uint8_t Leg_HangCalib(LegControl_t *leg);

/* 关节空间增量: 直接在模型角上给增量, 不经过逆解(起步阶段逆解在满伸处奇异, 只能用这个)。
 * dq1 / dq2 是模型角增量(rad, 腿平面内逆时针为正), 内部按 q1_direction / q2_direction 换算成电机角。
 * 例: 收腿 30 度 -> Leg_MoveJoint(&leg[0], +0.5236f, -0.5236f);
 *     镜像的那两条腿符号相反 -> Leg_MoveJoint(&leg[1], -0.5236f, +0.5236f);
 *     两种在实物上都表现为"每条腿两个电机的 theta 同向增大" 30 度。 */
uint8_t Leg_MoveJoint(LegControl_t *leg, float dq1, float dq2);

/* 读四条腿的反馈, 返回在线的腿数(0..4)。 */
uint8_t Leg_ReadAll(void);

#endif
