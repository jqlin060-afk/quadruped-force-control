#ifndef MOTOR_H
#define MOTOR_H

#define Pai 3.14159265358979323846f

#include <stdint.h>
#include "main.h"
#include "stm32h7xx_hal.h"
#include "usart.h"
#pragma pack(1)
typedef struct{
    uint8_t head[2];
    uint8_t idmode;
    int16_t tau;
	int16_t omega;
	int32_t  theta;
    uint16_t kp;
    uint16_t kd;
	uint16_t crc;
}MotorcmdRaw_t;

typedef struct {
    uint8_t id;
    float tau_nm;
    float omega_rad_s;
    float theta_rad;
    float kp;
    float kd;
} MotorcmdReal_t;

typedef struct{
	uint8_t head[2];
	uint8_t idmode;
    int16_t tau;
	int16_t omega;
	int32_t  theta;
	int8_t temp;
	uint16_t other;
	uint16_t crc;
}MotorbackRaw_t;

typedef struct {
    float tau;
    float omega;
    float theta;
    int8_t temp;
    uint8_t id;
    uint8_t status;
    /* 下面两个原来被丢在 MotorbackRaw_t::other 里，见电机手册表 2 第 96~110 位:
     *   bit 96-98   3bit  MERROR  0 正常 / 1 过热 / 2 过流 / 3 过压 / 4 编码器故障
     *   bit 99-110 12bit  FORCE   手册标"足端力 12bit 原始数据(0-4095)"，单位待确认
     *   bit 111      1bit  保留 */
    uint8_t  merror;
    uint16_t force_raw;
} MotorbackReal_t;
#pragma pack()

extern volatile MotorbackRaw_t motor_back;
extern volatile MotorbackReal_t Feedback_now;
extern volatile MotorbackReal_t Feedback[16];
extern volatile uint8_t Motor_Online[16];
/* 每台电机收到的回复计数。力控/调试时可以用它判断"这台还在不在回" */
extern volatile uint32_t Motor_RxCount[16];

uint16_t Motor_CRC16_CCITT(const uint8_t *data, uint16_t len);
void Motor_FillCmdCRC(MotorcmdRaw_t  *cmd);
uint8_t Motor_CheckBackCRC(MotorbackRaw_t *back);

void Motor_Receive_DMA(void);

void Motor_PollRx(void);

uint8_t Motor_IsOnline(uint8_t id);
void Motor_CheckOffline(void);

/* 单独问一台电机: 先把它标记成离线, 发一帧 0 力矩, 然后等它回一帧。
 * 半双工 485 上不能"连发多条命令再去收回复", 从站的回复会和主站正在发的
 * 下一条命令撞在一起, 结果就是有些电机永远读不到。所以必须一条一条问。
 * 返回 1 = 收到回复; 0 = 超时。 */
uint8_t Motor_Probe(uint8_t id, uint32_t timeout_ms);

MotorcmdRaw_t Motor_SetCmdReal(MotorcmdReal_t *Realcmd);
MotorbackReal_t Motor_ReceDataReal(MotorbackRaw_t *Rawback);
uint8_t Motor_CollectCmd(MotorcmdRaw_t *cmd);
HAL_StatusTypeDef Motor_SendCmdDMA(void);
#endif
