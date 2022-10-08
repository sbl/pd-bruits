#pragma once

#define br_minimum(x, y) (y < x ? y : x)
#define br_maximum(x, y) (x < y ? y : x)
#define br_clamp(x, minVal, maxVal) (br_minimum(br_maximum(x, minVal), maxVal))
