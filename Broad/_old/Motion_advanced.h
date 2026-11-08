#ifndef _MOTION_ADVANCED_H
#define _MOTION_ADVANCED_H

#include "Motor_motion.h"

uint8_t Leg_SetSwing(LegControl_t *leg,
                     float center_x,
                     float center_y,
                     float amplitude_x,
                     float amplitude_y,
                     float frequency_hz,
                     float phase_offset);

#endif