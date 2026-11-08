#ifndef LEG_FORCE_H
#define LEG_FORCE_H

#include "Motor_motion.h"

/* ==========================================================================
 * VMC 力控层（虚拟腿 + 纯力矩）
 *
 * 电机侧: kp 恒为 0, 只下发前馈力矩和速度阻尼 -> 电机是一个纯力矩源。
 * MCU 侧: 用"虚拟腿"算力矩（沿腿长 + 摆角各一条弹簧）:
 *
 *     沿腿长: F_r = kr   * (r_d - r)    + 重力分量
 *     摆角  : Q_b = kphi * (beta_d - beta) + 重力分量
 *
 *     alpha = (q1 - q2)/2,   Q_a = F_r * dr/dalpha
 *     tau1  = (Q_a + Q_b)/2,  tau2 = (Q_a - Q_b)/2
 *
 * 阻尼交给驱动器的 kd —— 它在自己的 20kHz 环里算, 没有通信延迟;
 * MCU 这边发命令到收反馈要 1.2~2ms, 用这么慢的速度量做微分会振。
 *
 * 同时这一层还负责观测: 由反馈 tau 反解足端力 (F = (J^T)^-1 * tau)。
 * ========================================================================== */

#define LF_G             9.81f
#define LF_MASS_DEFAULT  12.5f     /* 机体重量默认值(kg), 实测后覆盖 */
#define LF_NLEG          4u
#define LF_CM2M          0.01f     /* 运动学用 cm, 力和力矩换算用 m */

/* ===== 观测量（调试器 Watch 里直接加） ===== */
extern float   g_lf_fx[LF_NLEG];   /* 腿对地面的力: x 前后, N */
extern float   g_lf_fy[LF_NLEG];   /* 腿对地面的力: y 竖直, 承重时 < 0 */
extern float   g_lf_fz[LF_NLEG];   /* = -g_lf_fy, 承重时 > 0 */
extern float   g_lf_sum_fz;        /* 四腿竖直力之和; 站着应 ≈ m*g */
extern float   g_lf_mass_est;      /* = sum_fz / LF_G; 站着就是体重实测值 */
extern uint8_t g_lf_n_valid;       /* 上面两个数只统计有效的腿; < 4 时 mass_est 偏小 */
extern uint8_t g_lf_ok[LF_NLEG];   /* 0 = 该腿 J 奇异或离线, 力不可信 */
extern float   g_lf_jac_err;       /* 解析 J vs 数值 J 的最大相对误差 */
extern uint8_t g_lf_jac_bad;       /* 自检失败的位形个数 */

/* ===== 控制量 ===== */

/* 总开关。0 = 完全退回原来的位控(Leg_Update), 行为和改造前一样。 */
extern uint8_t g_lf_enable;

/* 只允许力控的腿: bit0=leg[0] ... bit3=leg[3]。第一次测试请用 0x01。 */
extern uint8_t g_lf_leg_mask;

/* 虚拟腿刚度 */
extern float   g_lf_kr;            /* 腿长刚度 N/m, 1000 = 1 N/mm */
extern float   g_lf_kphi;          /* 摆角刚度 N·m/rad */

/* 驱动器阻尼（转子侧）。**按相位分开**，因为两类腿的等效质量差 8 倍：
 *   摆动腿只有腿自己(~0.8kg)  -> 临界阻尼 kd≈0.027
 *   支撑腿扛着整机(~6.6kg)    -> 临界阻尼 kd≈0.075
 * 如果摆动腿抬不起来（滞后好厘米），就是这个值给太大了；落地弹就是太小。 */
extern float   g_lf_kd_swing;
extern float   g_lf_kd_stance;

/* 重力前馈缩放: 0 = 不发力(靠弹簧硬扛, 会下沉), 1.0 = 按体重理论值顶住。
 * 出厂 1.0: 体重由前馈承担, 弹簧只负责"软硬"和扰动。
 * 下沉量 = (单腿承重 − 前馈力) / kr。 */
extern float   g_lf_ff_scale;

/* 机体质量(kg)。0 = 用 LF_MASS_DEFAULT。 */
extern float   g_lf_mass;

/* tau 换算: 真实输出轴力矩 = 电机反馈值 × 这个数（已确认为 6.33 转子侧）。
 * 代码内部一律用输出侧 N·m 计算, 收发两端各自换算。 */
extern float   g_lf_tau_scale;

extern float   g_lf_tau_limit;     /* 单电机力矩上限, 输出侧 N·m */
extern float   g_lf_tau_step;      /* 每拍最大变化量 N·m */
extern float   g_lf_ff_max_frac;   /* 前馈单腿硬上限(占体重比例), 防 n_stance 算错 */

/* 触地判定阈值（用上一拍的足端竖直力 g_lf_fz 判断）。
 * 重力前馈只在**触地**时给：否则摆动末端会把腿压长 ff/kr（2 腿支撑时 6cm）才落地。 */
extern float   g_lf_contact_on;    /* 力超过这个值 -> 认为触地 */
extern float   g_lf_contact_off;   /* 力低于这个值 -> 认为离地（迟滞） */

/* 逐腿前馈微调（N），用来把四条腿的负载拉匀。
 * 软弹簧会把重心偏移放大成几厘米的高度差，靠这个给个别腿加减前馈把它抹平。
 * 正值 = 这条腿多出力。全 0 = 不启用。 */
extern float   g_lf_ff_trim[LF_NLEG];

/* 当前在支撑相的腿数(分摊体重用) */
extern float   g_lf_n_stance;

/* ===== 通信看门狗 =====
 * 电机掉线后不会变软: 它一直执行最后一帧, 变成一根刚性支柱。所以只要有一条腿
 * 掉线就强制退回位控。 */
extern uint8_t g_lf_comm_fault;    /* 1 = 有腿掉线(只读) */
extern uint8_t g_lf_comm_strict;   /* 1 = 掉线时自动把 g_lf_enable 清零(默认) */

/* 摆动相的逐腿刚度倍率。某条腿抬不起来（机械摩擦/装配紧）时，只把那条腿在
 * **摆动相**的 kr 放大，不影响它支撑时的承重。全 1.0 = 不启用；常用 1.5~2.5。 */
extern float   g_lf_swing_boost[LF_NLEG];

/* ===== 防翻转保护 =====
 * 五连杆在 r = 13.5cm(完全折叠) 和 r = 38.5cm(完全伸直) 两处退化。r -> 13.5 附近
 * det(J) -> 0, 腿撑不住载荷会被压过退化点, 之后脚端能沿半径 26cm 的整圆乱滑
 * —— 就是"脚从最低点直接翻到最高点卡死"。越界就把弹簧变硬, 把腿拖回来。 */
extern float   g_lf_r_min;         /* 安全下界 cm, 默认 17 */
extern float   g_lf_r_max;         /* 安全上界 cm, 默认 36 */
extern float   g_lf_err_limit;     /* 足端偏离目标上限 cm, 默认 8 */

/* 五连杆解析 Jacobian。
 * 输入 theta1/theta2 是电机角, 输出 J 行主序 (m/rad):
 *   J[0]=dx/dth1, J[1]=dx/dth2, J[2]=dy/dth1, J[3]=dy/dth2
 * point 可为 NULL。返回 0 = 位形奇异, J 不可用。 */
uint8_t Leg_Jacobian(const LegParam_t *param,
                     float theta1, float theta2,
                     float J[4], LegPoint_t *point);

/* 由反馈关节力矩反解足端力: J^T * F = tau  ==>  F = (J^T)^-1 * tau
 * (fx, fy) 是腿对地面的力, N。承重时 fy < 0。 */
uint8_t Leg_MeasuredFootForce(const LegControl_t *leg, float *fx, float *fy);

void Leg_ForceInit(void);          /* 上电一次: Jacobian 自检 */
void Leg_ForceObserveAll(void);    /* 每拍: 刷新观测量, 不控制 */

uint8_t Leg_ForceBeginCycle(void); /* 每拍先调: 统计支撑腿数 + 通信看门狗 */
uint8_t Leg_ForceUpdate(LegControl_t *leg);   /* VMC 版 Leg_Update */

#endif
