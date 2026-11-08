#include "Motor_advanced.h"

/* 位置平滑用的"上一次真正下发的角度"。上电时电机停在自己的零点(theta=0), 初值 0。 */
static float pos_ref[16] = {0.0f};

void Motor_SetZeroTorque(uint8_t id)
{
    MotorcmdReal_t real_cmd = {
        .id = id,
        .tau_nm = 0.0f,
        .omega_rad_s = 0.0f,
        .theta_rad = 0.0f,
        .kp = 0.0f,
        .kd = 0.0f
    };

    MotorcmdRaw_t raw_cmd = Motor_SetCmdReal(&real_cmd);
    Motor_CollectCmd(&raw_cmd);
}

void Motor_SetPosition(uint8_t id, float theta)
{
    float out;

    if (id >= 16) {
        return;
    }

    pos_ref[id] += (theta - pos_ref[id]) * pos_smooth;
    out = pos_ref[id];

    MotorcmdReal_t real_cmd = {
        .id = id,
        .tau_nm = 0.0f,
        .omega_rad_s = 0.0f,
        .theta_rad = out,
        .kp = 1.0f,
        .kd = 0.05f
    };

    MotorcmdRaw_t raw_cmd = Motor_SetCmdReal(&real_cmd);
    Motor_CollectCmd(&raw_cmd);
}

void Motor_SetTorque(uint8_t id, float tau, float kd)
{
    if (id >= 16) {
        return;
    }

    /* 手册量程 0 <= Kw <= 25.599, 超了就打满, 避免打包溢出 */
    if (kd < 0.0f)  { kd = 0.0f; }
    if (kd > 25.5f) { kd = 25.5f; }

    /* 手册力矩编码上限 |tau| <= 127.99 Nm(转子侧)。这里只做防溢出硬夹,
     * 真正的力矩限幅在 Leg_Force.c 里按 g_lf_tau_limit 做。 */
    if (tau > 120.0f)  { tau = 120.0f; }
    if (tau < -120.0f) { tau = -120.0f; }

    /* kp = 0, 所以 theta 不参与控制; 填当前实测角度, 万一以后 kp 不为 0
     * 也不会因为目标突变而跳。 */
    MotorcmdReal_t real_cmd = {
        .id = id,
        .tau_nm = tau,
        .omega_rad_s = 0.0f,
        .theta_rad = Feedback[id].theta,
        .kp = 0.0f,
        .kd = kd
    };

    /* 顺手把位置平滑器追到当前实际角度。这样运行中从力矩模式切回位置模式时,
     * 位置环不会从一个很旧的目标开始斜坡。
     * (g_lf_enable 是允许运行中改的, 所以这一句有用) */
    pos_ref[id] = Feedback[id].theta;

    MotorcmdRaw_t raw_cmd = Motor_SetCmdReal(&real_cmd);
    Motor_CollectCmd(&raw_cmd);
}

void Motor_SyncPosition(uint8_t id, float theta)
{
    if (id >= 16) {
        return;
    }

    pos_ref[id] = theta;
}

void Motor_ResetPositionControl(uint8_t id)
{
    if (id >= 16) {
        return;
    }

    pos_ref[id] = 0.0f;
}
