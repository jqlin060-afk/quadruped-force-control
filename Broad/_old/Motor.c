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
static HAL_StatusTypeDef Motor_TxStartNext(void);

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

    tx_busy = 0;
    Motor_TxStartNext();
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

static HAL_StatusTypeDef Motor_TxStartNext(void)
{
    if (tx_busy) {
        return HAL_BUSY;
    }

    for (uint8_t n = 0; n < 16U; n++) {
        uint8_t id = (uint8_t)((tx_next_id + n) % 16U);
        uint16_t bit = (uint16_t)(1U << id);

        if (tx_pending_mask & bit) {
            if (tx_buf[id].head[0] != 0xFE || tx_buf[id].head[1] != 0xEE) {
                tx_pending_mask &= (uint16_t)~bit;
                continue;
            }

            tx_active_cmd = tx_buf[id];
            tx_next_id = (uint8_t)((id + 1U) % 16U);
            tx_pending_mask &= (uint16_t)~bit;
            tx_busy = 1;

            HAL_StatusTypeDef status =
                HAL_UART_Transmit_DMA(&huart2,
                                      (uint8_t *)&tx_active_cmd,
                                      sizeof(tx_active_cmd));

            if (status != HAL_OK) {
                tx_busy = 0;
                tx_pending_mask |= bit;
            }

            return status;
        }
    }

    return HAL_OK;
}

HAL_StatusTypeDef Motor_SendCmdDMA(void)
{
    uint8_t  guard;
    uint32_t t0;

    for (guard = 0U; guard < 16U; guard++) {
        if (tx_pending_mask == 0U) {
            break;
        }

        if (tx_busy) {
            t0 = HAL_GetTick();
            while (tx_busy && (HAL_GetTick() - t0) < 2U) {
            }
            if (tx_busy) {
                tx_busy = 0;
            }
        }

        if (Motor_TxStartNext() != HAL_OK) {
            break;
        }
    }

    return HAL_OK;
}
