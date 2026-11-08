#ifndef MOTOR_ADVANCED_H
#define MOTOR_ADVANCED_H

#include "Motor.h"
#include <math.h>

#define accel_limit 2.0f

/* 位置指令一阶平滑系数: 每次调用朝目标挪这个比例。
 * 1.0 = 不平滑(直接用目标值); 0.10 左右比较温和。
 * 主循环 1ms 一次时, 0.10 约等于 10ms 时间常数。 */
#define pos_smooth 0.10f

void Motor_SetSpeed(uint8_t id, float omega, float alpha_now, float Km);
MotorbackReal_t Motor_GetFeedback(void);
void Motor_SetZeroTorque(uint8_t id);
void Motor_SetPosition(uint8_t id, float theta);
/* 力位混合下发（电机手册 §7.7）。驱动器实际执行:
 *     tau_out = tau + kp*(theta - theta_now) + kd*(0 - omega_now)
 * 传 kp=1.0 / kd=0.05 / tau=0 时与 Motor_SetPosition 完全等价。
 * kp、kd 会按手册量程 [0, 25.5] 夹紧。 */
void Motor_SetImpedance(uint8_t id, float theta, float tau, float kp, float kd);
/* 把位置平滑器对齐到当前真实角度(比如标定之后), 避免第一段从旧值斜坡过来。 */
void Motor_SyncPosition(uint8_t id, float theta);
void Motor_ResetPositionControl(uint8_t id);
void Motor_SetSpeedSmooth(uint8_t id, float omega_target, float dt);
float Motor_GetFilteredAlpha(uint8_t id, float dt);
void Motor_SetAngle(uint8_t id, float theta, float dt);
void Motor_AddAngle(uint8_t id, float delta_theta);
void Motor_SetSwing(uint8_t id,
                    float center,
                    float amplitude,
                    float frequency);
#endif
