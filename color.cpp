#include "color.h"

color_t colorFromUint32ARGB(uint32_t c) {
    return {
        static_cast<uint8_t>((c & 0x00FF0000) >> 16),
        static_cast<uint8_t>((c & 0x0000FF00) >> 8),
        static_cast<uint8_t>((c & 0x000000FF))
    };
}

color_t colorFromUint32RGBA(uint32_t c) {
    return {
        // static_cast<uint8_t>((c & 0xFF000000) >> 24),
        // static_cast<uint8_t>((c & 0x00FF0000) >> 16),
        // static_cast<uint8_t>((c & 0x0000FF00) >> 8)
        static_cast<uint8_t>((c & 0x000000FF)),
        static_cast<uint8_t>((c & 0x0000FF00) >> 8),
        static_cast<uint8_t>((c & 0x00FF0000) >> 16)
    };
}

color_t sumColors3(color_t c0, color_t c1, color_t c2) {
    return {
        static_cast<uint8_t>(clamp(c0.r + c1.r + c2.r, 255.0)),
        static_cast<uint8_t>(clamp(c0.g + c1.g + c2.g, 255.0)),
        static_cast<uint8_t>(clamp(c0.b + c1.b + c2.b, 255.0))
    };
}

color_t sumColors(color_t c0, color_t c1) {
    return {
        static_cast<uint8_t>(clamp(c0.r + c1.r, 255.0)),
        static_cast<uint8_t>(clamp(c0.g + c1.g, 255.0)),
        static_cast<uint8_t>(clamp(c0.b + c1.b, 255.0))
    };
}

color_t colorFromFloats(float r, float g, float b) {
    return {
        static_cast<uint8_t>(clamp(r * 255.0, 255.0)),
        static_cast<uint8_t>(clamp(g * 255.0, 255.0)),
        static_cast<uint8_t>(clamp(b * 255.0, 255.0))
    };
}