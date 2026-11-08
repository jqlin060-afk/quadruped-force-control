#include "Motor_advanced.h"

float omega_last[16] = {0.0f};
static float omega_ref[16] = {0.0f};
static float alpha_filt[16] = {0.0f};
static float software_theta[16] = {0.0f};
/* 位置平滑用: 上一次真正下发的角度。每来一个新目标, 只朝它挪 pos_smooth 的比例。
 * 上电时电机就停在自己零点(theta=0), 所以初值取 0 是对的。 */
static float pos_ref[16] = {0.0f};
static uint8_t swing_started[16] = {0};
static uint8_t swing_direction[16] = {0};
static uint32_t swing_next_tick[16] = {0};

void Motor_SetSpeed(uint8_t id, float omega, float alpha_now, float Km)
{
    MotorcmdReal_t real_cmd = {
        .id = id,
        .tau_nm = -Km * alpha_now,
        .omega_rad_s = omega,
        .theta_rad = 0.0f,
        .kp = 0.0f,
        .kd = 0.08f
    };

    MotorcmdRaw_t raw_cmd = Motor_SetCmdReal(&real_cmd);
    Motor_CollectCmd(&raw_cmd);
}

MotorbackReal_t Motor_GetFeedback(void)
{
    return Feedback_now;
}

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

    /* 一阶平滑: 不直接把目标丢给电机, 而是每次朝目标挪一点。
     * 好处: 目标突然跳很远(比如刚上电、或换目标点)时, 电机是斜坡过去而不是猛冲。
     * 代价: 有轻微滞后(主循环 1ms 一次时约 10ms), 对 1Hz 步态完全够用。 */
    pos_ref[id] += (theta - pos_ref[id]) * pos_smooth;
    out = pos_ref[id];

    software_theta[id] = out;

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

/* 力位混合下发。和 Motor_SetPosition 走同一套位置平滑器 pos_ref，
 * 区别只是 kp/kd 由调用者给、并且多一个前馈力矩 tau。 */
void Motor_SetImpedance(uint8_t id, float theta, float tau, float kp, float kd)
{
    float out;

    if (id >= 16) {
        return;
    }

    /* 手册量程 0 <= Kpos <= 25.599，超了就打满，避免打包溢出 */
    if (kp < 0.0f)      { kp = 0.0f; }
    if (kp > 25.5f)     { kp = 25.5f; }
    if (kd < 0.0f)      { kd = 0.0f; }
    if (kd > 25.5f)     { kd = 25.5f; }

    /* 手册力矩编码上限 |tau| <= 127.99 Nm，这里只做防溢出的硬夹，
     * 真正的力矩限幅在 Leg_Force.c 里按 g_lf_tau_limit 做。 */
    if (tau > 120.0f)   { tau = 120.0f; }
    if (tau < -120.0f)  { tau = -120.0f; }

    pos_ref[id] += (theta - pos_ref[id]) * pos_smooth;
    out = pos_ref[id];

    software_theta[id] = out;

    MotorcmdReal_t real_cmd = {
        .id = id,
        .tau_nm = tau,
        .omega_rad_s = 0.0f,
        .theta_rad = out,
        .kp = kp,
        .kd = kd
    };

    MotorcmdRaw_t raw_cmd = Motor_SetCmdReal(&real_cmd);
    Motor_CollectCmd(&raw_cmd);
}

void Motor_SyncPosition(uint8_t id, float theta)
{
    if (id >= 16) {
        return;
    }

    pos_ref[id] = theta;
    software_theta[id] = theta;
}
void Motor_ResetPositionControl(uint8_t id)
{
    if (id >= 16) {
        return;
    }

    omega_ref[id] = 0.0f;
    omega_last[id] = 0.0f;
    alpha_filt[id] = 0.0f;
    software_theta[id] = 0.0f;
    pos_ref[id] = 0.0f;
    swing_started[id] = 0;
    swing_direction[id] = 0;
    swing_next_tick[id] = 0;
}

void Motor_SetSpeedSmooth(uint8_t id, float omega_target, float dt)
{
    float max_delta;
    float delta;
    float alpha_now;

    if (id >= 16 || dt <= 0.0f) {
        return;
    }

    max_delta = accel_limit * dt;
    delta = omega_target - omega_ref[id];

    if (delta > max_delta) {
        delta = max_delta;
    }

    if (delta < -max_delta) {
        delta = -max_delta;
    }

    alpha_now = Motor_GetFilteredAlpha(id, dt);
    omega_ref[id] += delta;

    Motor_SetSpeed(id, omega_ref[id], alpha_now, 0.1f);
}

float Motor_GetFilteredAlpha(uint8_t id, float dt)
{
    float omega_now;
    float alpha_raw;
    float k = 0.1f;

    if (id >= 16 || dt <= 0.0f) {
        return 0.0f;
    }

    omega_now = Feedback[id].omega;
    alpha_raw = (omega_now - omega_last[id]) / dt;
    omega_last[id] = omega_now;

    alpha_filt[id] = alpha_filt[id] + k * (alpha_raw - alpha_filt[id]);

    return alpha_filt[id];
}

void Motor_SetAngle(uint8_t id, float theta, float dt)
{
    (void)dt;

    if (id >= 16) {
        return;
    }

    swing_started[id] = 0;
    Motor_SetPosition(id, theta);
}

void Motor_AddAngle(uint8_t id, float delta_theta)
{
    if (id >= 16) {
        return;
    }

    swing_started[id] = 0;
    software_theta[id] += delta_theta;
    Motor_SetPosition(id, software_theta[id]);
}

void Motor_SetSwing(uint8_t id,
                    float center,
                    float amplitude,
                    float frequency)
{
    uint32_t now;
    uint32_t half_period_ms;

    if (id >= 16 || frequency <= 0.0f) {
        return;
    }

    if (amplitude < 0.0f) {
        amplitude = -amplitude;
    }

    half_period_ms = (uint32_t)(500.0f / frequency);
    if (half_period_ms == 0U) {
        half_period_ms = 1U;
    }

    now = HAL_GetTick();

    if (!swing_started[id]) {
        swing_started[id] = 1;
        swing_direction[id] = 1;
        swing_next_tick[id] = now + half_period_ms;
    }

    while ((int32_t)(now - swing_next_tick[id]) >= 0) {
        swing_direction[id] = !swing_direction[id];
        swing_next_tick[id] += half_period_ms;
    }

    if (swing_direction[id]) {
        software_theta[id] = center + amplitude;
    } else {
        software_theta[id] = center - amplitude;
    }

    Motor_SetPosition(id, software_theta[id]);
}
