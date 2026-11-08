/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2024 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "dma.h"
#include "usart.h"
#include "gpio.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include "Motor.h"
#include "Motor_advanced.h"
#include "Motor_motion.h"
#include "Motion_advanced.h"      /* 运动层: 步态 */
#include "Leg_Force.h"            /* 力控层: Jacobian / 足端力观测（目前只观测，不控制） */
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
/* motor1 必须是"左电机", motor2 必须是"右电机" */

/* ===== 走路参数 ===== */
#define STEP_AX      6.0f     /* 足端前后摆幅 (cm)。整机速度 ≈ 2*AX/(0.6*T) ≈ 10.7 cm/s */
#define STEP_AY      5.0f     /* 抬腿高度 (cm) */
#define STEP_FREQ    0.8f     /* 步频 (Hz) */

#define WALK_EN       1      /* 1 = 收腿后直接走; 0 = 只站立 */

/* 力控总开关(编译期)。1 = 把柔顺控制的代码路径编进来，但**是否真的生效**由运行期的
 * g_lf_enable 决定（默认 0）。也就是说: 烧进去的行为和现在完全一样，
 * 想看柔顺效果就在调试器里把 g_lf_enable 改成 1，不用重编译。
 * 0 = 连代码都不编进来。 */
#define FORCE_EN     1
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */

/* 调试器里看这四个: 1 = 该腿标定成功(会上报、会动), 0 = 失败(被禁用, 不收位置命令) */
volatile uint8_t Leg_CalibOK[4];
/* 握手阶段实际读到的在线腿数 + 握手扫了多少遍, 用来定位是哪一步出问题 */
volatile uint8_t  Leg_OnlineNum;
volatile uint16_t Leg_HandshakeTry;
/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MPU_Config(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{

  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */
  MPU_Config();
  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_USART2_UART_Init();
  /* USER CODE BEGIN 2 */
  Motor_Receive_DMA();

  /* 等电机启动完, 否则最早几条命令会被丢掉 */
  HAL_Delay(1000);

  /* ============ 悬垂上电标定 + 缓慢收腿 ============
   * 前提: 狗悬空、四条腿自然垂直朝下、静止不动。 */

  /* 1) 建四条腿的模型 */
  Led_ControlUnitInit();

  /* 1b) 力控自检: 解析 Jacobian 对数值差分逐点校验，结果看 g_lf_jac_err / g_lf_jac_bad。
   *     正常应该是 g_lf_jac_err < 1e-2 且 g_lf_jac_bad == 0。 */
  Leg_ForceInit();

  /* 2) 握手: 电机只有收到命令才会回一帧反馈, 光轮询是读不到 theta 的。
   *    这里发 0 力矩去"问"一帧: 不发力, 腿保持自然悬垂, 正好是标定要的姿态。
   *    半双工 485 上不能连发多条命令再去收回复(从站的回复会和主站正在发的
   *    下一条命令撞车), 所以必须一条一条问: 发一帧 -> 等它回 -> 再问下一个。
   *    每个电机最多扫 3 遍。 */
  {
    uint8_t sweep;

    Leg_OnlineNum = 0U;

    for (sweep = 0U; sweep < 3U; sweep++) {
      uint8_t i;

      Leg_HandshakeTry = sweep;

      for (i = 0U; i < 4U; i++) {
        uint8_t ok1 = Motor_IsOnline(leg[i].param.motor1_id);
        uint8_t ok2 = Motor_IsOnline(leg[i].param.motor2_id);

        if (!ok1) {
          Motor_Probe(leg[i].param.motor1_id, 6U);
        }
        if (!ok2) {
          Motor_Probe(leg[i].param.motor2_id, 6U);
        }
      }

      Leg_OnlineNum = Leg_ReadAll();
      if (Leg_OnlineNum >= 4U) {
        break;
      }
    }

    /* 让腿在这段时间里松着、自然垂到底 */
    Motor_SetZeroTorque(leg[0].param.motor1_id);
    Motor_SetZeroTorque(leg[0].param.motor2_id);
    Motor_SetZeroTorque(leg[1].param.motor1_id);
    Motor_SetZeroTorque(leg[1].param.motor2_id);
    Motor_SetZeroTorque(leg[2].param.motor1_id);
    Motor_SetZeroTorque(leg[2].param.motor2_id);
    Motor_SetZeroTorque(leg[3].param.motor1_id);
    Motor_SetZeroTorque(leg[3].param.motor2_id);
    Motor_SendCmdDMA();
    HAL_Delay(300U);
    Motor_PollRx();
  }

  /* 3) 悬垂标定: 把当前姿态声明成"垂直朝下、满伸", 定零位,
   *    并把下发值改成绝对角(q_off = 当前 theta)。
   *    标定失败的那条腿直接禁用(ready = 0), 它的 Leg_MoveTo/Leg_Update 会
   *    立刻返回 0、不发任何位置命令 —— 也就是说那条腿保持松着不动, 而不是拿
   *    没标定的模型骗它。其余能标定的腿照常工作, 不会因为一条腿拖死整机。 */
  {
    uint8_t i;

    for (i = 0U; i < 4U; i++) {
      if (Leg_HangCalib(&leg[i])) {
        Leg_CalibOK[i] = 1U;
      } else {
        Leg_CalibOK[i] = 0U;
        leg[i].ready = 0U;          /* 禁用这条腿 */
      }
    }
  }

  /* 3) 试验 Leg_MoveTo: 足端目标从满伸(MMA_HANG_Y = -38.5)平滑抬到 -26。
   * 全程走运动学逆解, 用来验证:
   *   四条腿是否同步收上来、足端是否只在轴线上下动(不往一侧偏)、
   *   收完四条腿高度是否一致(看 leg[i].real_point)。
   * 分 400 步是为了不让电机一步跳 76 度(满伸 -> r=26 需要 alpha 76.1 度)。 */
  {
    const uint32_t   STEPS = 400U;
    const float      TARGET_Y = -26.0f;
    uint32_t n;

    for (n = 0U; n <= STEPS; n++) {
      uint8_t i;
      float   t = (float)n / (float)STEPS;
      float   y = MMA_HANG_Y + (TARGET_Y - MMA_HANG_Y) * t;

      for (i = 0U; i < 4U; i++) {
        Leg_MoveTo(&leg[i], 0.0f, y);
      }

      Motor_PollRx();
      Motor_CheckOffline();
      HAL_Delay(3);
    }
  }

  /* 4) 收腿之后: 按 WALK_EN 选"站立"还是"走路" */
#if WALK_EN
  /* trot: 对角腿同相。左前+右后 一组(相位 0), 右前+左后 一组(相位 Pai)。
   * 哪条腿前后方向反了, 就把那条腿的 STEP_AX 改成负的。 */
  Leg_SetSwing(&leg[0], 0.0f, STAND_Y, STEP_AX, STEP_AY, STEP_FREQ, 0.0f);
  Leg_SetSwing(&leg[3], 0.0f, STAND_Y, STEP_AX, STEP_AY, STEP_FREQ, 0.0f);
  Leg_SetSwing(&leg[1], 0.0f, STAND_Y, STEP_AX, STEP_AY, STEP_FREQ, Pai);
  Leg_SetSwing(&leg[2], 0.0f, STAND_Y, STEP_AX, STEP_AY, STEP_FREQ, Pai);
#else
  /* 只站立: 四条腿保持在 (0, STAND_Y) */
  {
    uint8_t i;

    for (i = 0U; i < 4U; i++) {
      leg[i].hold_x = 0.0f;
      leg[i].hold_y = STAND_Y;
    }
  }
#endif
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
#if FORCE_EN
	  /* 先统计这一拍有几条腿在支撑相(体重按这个数分摊) */
	  Leg_ForceBeginCycle();

	  Leg_ForceUpdate(&leg[0]);
	  Leg_ForceUpdate(&leg[1]);
	  Leg_ForceUpdate(&leg[2]);
	  Leg_ForceUpdate(&leg[3]);
#else
	  Leg_Update(&leg[0]);
	  Leg_Update(&leg[1]);
	  Leg_Update(&leg[2]);
	  Leg_Update(&leg[3]);
#endif
	  
//	  Motor_SetZeroTorque(2);
//      Motor_SendCmdDMA();
//	  Motor_SetZeroTorque(3);
//      Motor_SendCmdDMA();
//	  Motor_SetZeroTorque(4);
//      Motor_SendCmdDMA();
//	  Motor_SetZeroTorque(5);
//      Motor_SendCmdDMA();
//	  Motor_SetZeroTorque(6);
//      Motor_SendCmdDMA();
//	  Motor_SetZeroTorque(7);
//      Motor_SendCmdDMA();
//	  Motor_SetZeroTorque(8);
//      Motor_SendCmdDMA();
//	  Motor_SetZeroTorque(9);
      Motor_SendCmdDMA();

	  Motor_PollRx();
	  Motor_CheckOffline();

	  Leg_ReadFeedback(&leg[0]);
	  Leg_ReadFeedback(&leg[1]);
	  Leg_ReadFeedback(&leg[2]);
	  Leg_ReadFeedback(&leg[3]);

	  /* 只观测, 不控制: 由 tau 反解足端力, 刷新 g_lf_fz[] / g_lf_sum_fz / g_lf_mass_est */
	  Leg_ForceObserveAll();

	  HAL_Delay(1);
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Supply configuration update enable
  */
  HAL_PWREx_ConfigSupply(PWR_LDO_SUPPLY);

  /** Configure the main internal regulator output voltage
  */
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE3);

  while(!__HAL_PWR_GET_FLAG(PWR_FLAG_VOSRDY)) {}

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_DIV1;
  RCC_OscInitStruct.HSICalibrationValue = 64;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_NONE;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2
                              |RCC_CLOCKTYPE_D3PCLK1|RCC_CLOCKTYPE_D1PCLK1;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_HSI;
  RCC_ClkInitStruct.SYSCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_HCLK_DIV1;
  RCC_ClkInitStruct.APB3CLKDivider = RCC_APB3_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_APB1_DIV1;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_APB2_DIV1;
  RCC_ClkInitStruct.APB4CLKDivider = RCC_APB4_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_1) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE BEGIN 4 */

/* USER CODE END 4 */

/**
  * @brief  Configure the MPU attributes
  * @retval None
  */
static void MPU_Config(void)
{
  MPU_Region_InitTypeDef MPU_InitStruct = {0};

  /* Disables the MPU */
  HAL_MPU_Disable();

  /** Initializes and configures the Region and the memory to be protected
  */
  MPU_InitStruct.Enable = MPU_REGION_ENABLE;
  MPU_InitStruct.Number = MPU_REGION_NUMBER0;
  MPU_InitStruct.BaseAddress = 0x0;
  MPU_InitStruct.Size = MPU_REGION_SIZE_4GB;
  MPU_InitStruct.SubRegionDisable = 0x87;
  MPU_InitStruct.TypeExtField = MPU_TEX_LEVEL0;
  MPU_InitStruct.AccessPermission = MPU_REGION_NO_ACCESS;
  MPU_InitStruct.DisableExec = MPU_INSTRUCTION_ACCESS_DISABLE;
  MPU_InitStruct.IsShareable = MPU_ACCESS_SHAREABLE;
  MPU_InitStruct.IsCacheable = MPU_ACCESS_NOT_CACHEABLE;
  MPU_InitStruct.IsBufferable = MPU_ACCESS_NOT_BUFFERABLE;

  HAL_MPU_ConfigRegion(&MPU_InitStruct);

  /* Enables the MPU */
  HAL_MPU_Enable(MPU_PRIVILEGED_DEFAULT);

}

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
