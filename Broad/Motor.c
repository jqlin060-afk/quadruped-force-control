#include "Motor.h"
#include <string.h>

extern UART_HandleTypeDef huart2;

volatile MotorbackRaw_t motor_back;
volatile MotorbackReal_t Feedback_now;
volatile MotorbackReal_t Feedback[16];
volatile uint8_t Motor_Online[16];
static uint32_t Motor_LastRxTick[16];

static uint8_t  rx_buf[256];
static uint16_t rx_read_pos = 0;

static uint8_t  rx_frame[16];

static uint8_t  rx_len = 0;

static MotorcmdRaw_t tx_buf[16];
static MotorcmdRaw_t tx_active_cmd;
static volatile uint16_t tx_pending_mask = 0;
static volatile uint8_t  tx_next_id = 0;
static volatile uint8_t  tx_busy = 0;

/* 每台电机收到的回复计数。TDMA 用它判断"这一台的回复到了没有" */
volatile uint32_t Motor_RxCount[16];

/* 一帧命令发出后等这台电机回复的时间上限(微秒)。
 * 太小会误判超时、太大主循环变慢。先给 300。 */
#define MOTOR_SLOT_US   300U

/* DWT 周期计数器, 做微秒级超时用(HAL_GetTick 只有 1ms 分辨率, 不够) */
static uint8_t  dwt_ready = 0;

static void Motor_InitCycleCounter(void)
{
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    DWT->CYCCNT = 0;
    DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
    dwt_ready = 1;
}

static uint32_t Motor_Cyc(void)
{
    return DWT->CYCCNT;
}

uint16_t Motor_CRC16_CCITT(const uint8_t *data, uint16_t len)
{
    uint16_t crc = 0x0000;

    while (len--) {
        crc ^= *data++;
        for (uint8_t i = 0; i < 8; i++) {
            if (crc & 0x0001) {
                crc = (crc >> 1) ^ 0x8408;
            } else {
                crc >>= 1;
            }
        }
    }

    return crc;
}

void Motor_FillCmdCRC(MotorcmdRaw_t *cmd)
{
    cmd->crc = Motor_CRC16_CCITT((uint8_t *)cmd, sizeof(MotorcmdRaw_t) - 2);
}

uint8_t Motor_CheckBackCRC(MotorbackRaw_t *back)
{
    uint16_t crc = Motor_CRC16_CCITT((uint8_t *)back, sizeof(MotorbackRaw_t) - 2);
    return crc == back->crc;
}

MotorbackReal_t Motor_ReceDataReal(MotorbackRaw_t *Rawback)
{
    MotorbackReal_t Realback;
    Realback.id = Rawback->idmode & 0x0F;
    Realback.status = (Rawback->idmode >> 4) & 0x07;
    Realback.tau = (float)Rawback->tau / 256.0f;
    Realback.omega = (float)Rawback->omega / 256.0f * 6.2831853f / 6.33f;
    Realback.theta = (float)Rawback->theta / 32768.0f * 6.2831853f / 6.33f;
    Realback.temp = Rawback->temp;
    /* other 是小端读出来的 16bit，对应帧内 bit 96~111 */
    Realback.merror = (uint8_t)(Rawback->other & 0x0007U);
    Realback.force_raw = (uint16_t)((Rawback->other >> 3) & 0x0FFFU);
    return Realback;
}

MotorcmdRaw_t Motor_SetCmdReal(MotorcmdReal_t *Realcmd)
{
    MotorcmdRaw_t Rawcmd;
    Rawcmd.head[0] = 0xFE;
    Rawcmd.head[1] = 0xEE;
    Rawcmd.idmode = (uint8_t)((Realcmd->id & 0x0F) | (1 << 4));
    Rawcmd.tau = (int16_t)(Realcmd->tau_nm * 256.0f);
    Rawcmd.omega = (int16_t)(Realcmd->omega_rad_s * 6.33f / 6.2831853f * 256.0f);
    Rawcmd.theta = (int32_t)(Realcmd->theta_rad * 6.33f / 6.2831853f * 32768.0f);
    Rawcmd.kp = (uint16_t)(Realcmd->kp * 1280.0f);
    Rawcmd.kd = (uint16_t)(Realcmd->kd * 1280.0f);
    Motor_FillCmdCRC(&Rawcmd);
    return Rawcmd;
}

static void Motor_SaveFeedback(uint8_t *frame)
{
    MotorbackRaw_t back;

    memcpy(&back, frame, sizeof(back));

    if (back.head[0] != 0xFD || back.head[1] != 0xEE) {
        return;
    }

    if (!Motor_CheckBackCRC(&back)) {
        return;
    }

    MotorbackReal_t real = Motor_ReceDataReal(&back);

    if (real.id > 14U) {

        return;
    }

    motor_back = back;
    Feedback_now = real;
    Feedback[real.id] = real;
    Motor_LastRxTick[real.id] = HAL_GetTick();
    Motor_Online[real.id] = 1;
    Motor_RxCount[real.id]++;
}

static void Motor_RxByte(uint8_t byte)
{
    if (rx_len == 0U) {
        if (byte == 0xFD) {
            rx_frame[0] = byte;
            rx_len = 1U;
        }
        return;
    }

    if (rx_len == 1U) {
        if (byte == 0xEE) {
            rx_frame[1] = byte;
            rx_len = 2U;
        } else if (byte == 0xFD) {
            rx_frame[0] = byte;

        } else {
            rx_len = 0U;
        }
        return;
    }

    rx_frame[rx_len++] = byte;

    if (rx_len >= 16U) {
        Motor_SaveFeedback(rx_frame);
        rx_len = 0U;
    }
}

void Motor_Receive_DMA(void)
{
    rx_read_pos = 0U;
    rx_len = 0U;

    HAL_UARTEx_ReceiveToIdle_DMA(&huart2, rx_buf, sizeof(rx_buf));
}

void Motor_PollRx(void)
{
    uint16_t dma_pos;

    if (huart2.hdmarx == NULL) {
        return;
    }

    dma_pos = (uint16_t)(sizeof(rx_buf) - __HAL_DMA_GET_COUNTER(huart2.hdmarx));
    dma_pos %= sizeof(rx_buf);

    while (rx_read_pos != dma_pos) {
        Motor_RxByte(rx_buf[rx_read_pos]);
        rx_read_pos++;
        if (rx_read_pos >= sizeof(rx_buf)) {
            rx_read_pos = 0U;
        }
    }
}

uint8_t Motor_IsOnline(uint8_t id)
{
    if (id >= 16U) {
        return 0U;
    }
    return Motor_Online[id];
}

void Motor_CheckOffline(void)
{
    uint32_t now = HAL_GetTick();

    for (uint8_t i = 0; i < 16U; i++) {
        if (Motor_LastRxTick[i] != 0U &&
            (now - Motor_LastRxTick[i]) > 100U) {
            Motor_Online[i] = 0U;
        }
    }
}

uint8_t Motor_Probe(uint8_t id, uint32_t timeout_ms)
{
    MotorcmdReal_t  real_cmd;
    MotorcmdRaw_t   raw_cmd;
    uint32_t        t0;

    if (id >= 16U) {
        return 0U;
    }

    Motor_Online[id] = 0U;          /* 先清掉旧状态, 只看这一次有没有回 */

    real_cmd.id = id;
    real_cmd.tau_nm = 0.0f;
    real_cmd.omega_rad_s = 0.0f;
    real_cmd.theta_rad = 0.0f;
    real_cmd.kp = 0.0f;
    real_cmd.kd = 0.0f;

    raw_cmd = Motor_SetCmdReal(&real_cmd);
    Motor_CollectCmd(&raw_cmd);
    Motor_SendCmdDMA();

    t0 = HAL_GetTick();

    while ((HAL_GetTick() - t0) < timeout_ms) {
        Motor_PollRx();
        if (Motor_Online[id] != 0U) {
            return 1U;
        }
    }

    return 0U;
}

void HAL_UART_ErrorCallback(UART_HandleTypeDef *huart)
{
    if (huart->Instance != USART2) {
        return;
    }

    __HAL_UART_CLEAR_FLAG(huart, UART_CLEAR_OREF | UART_CLEAR_FEF |
                                UART_CLEAR_NEF | UART_CLEAR_PEF);

    huart->RxState = HAL_UART_STATE_READY;
    huart->ReceptionType = HAL_UART_RECEPTION_STANDARD;

    rx_len = 0U;
    rx_read_pos = 0U;

    HAL_UARTEx_ReceiveToIdle_DMA(huart, rx_buf, sizeof(rx_buf));
}

void HAL_UART_TxCpltCallback(UART_HandleTypeDef *huart)
{
    if (huart->Instance != USART2) {
        return;
    }

    /* 只清标志, **不要**在这里启动下一帧。
     * 下一帧什么时候发由 Motor_SendCmdDMA 的 TDMA 节奏决定(要等这一台的回复)。 */
    tx_busy = 0;
}

uint8_t Motor_CollectCmd(MotorcmdRaw_t *cmd)
{
    uint8_t id;

    if (cmd == NULL) {
        return 0;
    }

    id = cmd->idmode & 0x0F;

    uint32_t primask = __get_PRIMASK();

    __disable_irq();
    tx_buf[id] = *cmd;
    tx_pending_mask |= (uint16_t)(1U << id);

    if (!primask) {
        __enable_irq();
    }

    return 1;
}

/* ==========================================================================
 * 严格 TDMA: 发一帧 -> 等这台电机的回复(或超时) -> 再发下一帧。
 *
 * 为什么必须这样: 485 是半双工, 发送时 DE 把双绞线占住, 而电机同时在回复也在
 * 驱动同一条线 —— 两边电平一起被搞坏。于是:
 *   1) 电机收到的那帧 CRC 校验失败 -> 直接忽略 -> 继续执行上一条命令
 *      (那条是 kp=1.0 的位控帧) -> 表现就是"这个电机还是硬的";
 *   2) MCU 也收不到回复 -> 那条腿的反馈是陈旧的。
 * 8 帧连发时大约一半会撞掉, 所以看起来像"偶数硬、奇数软"。
 * ========================================================================== */
HAL_StatusTypeDef Motor_SendCmdDMA(void)
{
    uint8_t  n;
    uint16_t k;
    uint32_t slot_cyc;
    uint32_t t0;

    if (!dwt_ready) {
        Motor_InitCycleCounter();
    }

    /* 一台电机的槽时间: 发 42.5us + 回复延迟 + 收 40us */
    slot_cyc = (SystemCoreClock / 1000000U) * MOTOR_SLOT_US;
    if (slot_cyc == 0U) {
        slot_cyc = 64U * MOTOR_SLOT_US;
    }

    for (n = 0U; n < 16U; n++) {
        uint8_t  id = 0xFFU;
        uint16_t bit;
        uint32_t seq;

        if (tx_pending_mask == 0U) {
            break;
        }

        /* 从 tx_next_id 起轮询下一个待发的电机 */
        for (k = 0U; k < 16U; k++) {
            uint8_t cand = (uint8_t)((tx_next_id + k) % 16U);
            if (tx_pending_mask & (uint16_t)(1U << cand)) {
                id = cand;
                break;
            }
        }
        if (id == 0xFFU) {
            break;
        }

        bit = (uint16_t)(1U << id);

        if (tx_buf[id].head[0] != 0xFE || tx_buf[id].head[1] != 0xEE) {
            tx_pending_mask &= (uint16_t)~bit;
            continue;
        }

        /* 等上一帧真的发完 */
        t0 = Motor_Cyc();
        while (tx_busy && ((Motor_Cyc() - t0) < slot_cyc)) {
        }
        tx_busy = 0;

        tx_active_cmd = tx_buf[id];
        tx_next_id = (uint8_t)((id + 1U) % 16U);
        tx_pending_mask &= (uint16_t)~bit;

        seq = Motor_RxCount[id];        /* 记住发送前的回复计数 */
        tx_busy = 1;

        if (HAL_UART_Transmit_DMA(&huart2, (uint8_t *)&tx_active_cmd,
                                  sizeof(tx_active_cmd)) != HAL_OK) {
            tx_busy = 0;
            continue;
        }

        /* 等这台电机的回复。收到了马上发下一帧, 没收到就等满一个槽 */
        t0 = Motor_Cyc();
        while ((Motor_Cyc() - t0) < slot_cyc) {
            Motor_PollRx();
            if (Motor_RxCount[id] != seq) {
                break;
            }
        }
    }

    return HAL_OK;
}
