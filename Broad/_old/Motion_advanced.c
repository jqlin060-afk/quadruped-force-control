#include "Motion_advanced.h"

uint8_t Leg_SetSwing(LegControl_t *leg,
                     float center_x,
                     float center_y,
                     float amplitude_x,
                     float amplitude_y,
                     float frequency_hz,
                     float phase_offset)
{
    if (leg == 0 || !leg->ready || frequency_hz <= 0.0f) {
        return 0;
    }

    leg->mode = LEG_MODE_SWING;
    leg->swing_cx = center_x;
    leg->swing_cy = center_y;
    leg->swing_ax = amplitude_x;
    leg->swing_ay = amplitude_y;
    leg->swing_freq = frequency_hz;
    leg->swing_phase = phase_offset;
    leg->swing_t0 = HAL_GetTick();

    return 1;
}

