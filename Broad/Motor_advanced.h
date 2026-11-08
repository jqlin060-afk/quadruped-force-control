#ifndef MOTOR_ADVANCED_H
#define MOTOR_ADVANCED_H

#include "Motor.h"
#include <math.h>

/* 位置指令一阶平滑系数: 每次调用朝目标挪这个比例。
 * 1.0 = 不平滑(直接用目标值); 0.10 左右比较温和。
 * 主循环 1ms 一次时, 0.10 约等于 10ms 时间常数。 */
#define pos_smooth 0.10f

/* 零力矩: 电机主动产生力矩抵消自身摩擦(手册 §7.6), 腿保持"松"的状态。
 * 握手探测和上电悬垂阶段用。 */
void Motor_SetZeroTorque(uint8_t id);

/* 位置控制: 驱动器 kp=1.0 / kd=0.05 的位置环, 位控模式用。 */
void Motor_SetPosition(uint8_t id, float theta);

/* 纯力矩: kp 恒为 0, 只下发前馈力矩 + 速度阻尼。
 * 电机变成一个纯力矩源, 刚度全部由上层(Leg_Force 的虚拟腿)自己算。
 * 手册 §7.4/§7.5: 实际执行 tau_out = tau + kd*(0 - omega)。
 * kd 由驱动器在 20kHz 环里算, 没有通信延迟 —— 所以阻尼放这里, 位置反馈放上层。 */
void Motor_SetTorque(uint8_t id, float tau, float kd);

/* 把位置平滑器对齐到当前真实角度(标定之后用), 避免第一段从旧值斜坡过来。 */
void Motor_SyncPosition(uint8_t id, float theta);

void Motor_ResetPositionControl(uint8_t id);

#endif
