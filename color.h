#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>

typedef struct color_t {
  uint8_t r;
  uint8_t g;
  uint8_t b;
} color_t;

const color_t COLOR_WHITE      = {255, 255, 255};
const color_t COLOR_BLACK      = {0,   0,   0  };
const color_t COLOR_GREEN      = {0,   255, 0  };
const color_t COLOR_BLUE       = {0,   0,   255};
const color_t COLOR_RED        = {255, 0,   0  };
const color_t COLOR_YELLOW     = {255, 255, 0  };
const color_t COLOR_PURPLE     = {255, 0,   255};
const color_t COLOR_CYAN       = {0,   255, 255};

inline uint32_t colorToUint32(color_t c) {
    return 0x00000000 | (c.r << 16) | (c.g << 8) |  c.b;    
}

color_t colorFromUint32ARGB(uint32_t c);
color_t colorFromUint32RGBA(uint32_t c);
inline float clamp(float v, float max) {
    return v > max ? max : v;
}

inline color_t mulScalarColor(double x, color_t color) {
    color_t result = {
        static_cast<uint8_t>(clamp(x * color.r, 255.0)),
        static_cast<uint8_t>(clamp(x * color.g, 255.0)),
        static_cast<uint8_t>(clamp(x * color.b, 255.0))
    };

    return result;
}

color_t sumColors3(color_t c0, color_t c1, color_t c2);
color_t sumColors(color_t c0, color_t c1);
color_t colorFromFloats(float r, float g, float b);

#endif