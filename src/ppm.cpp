#include "ppm.h"

const double antiAliasMatrix[SAMPLE_RATE][SAMPLE_RATE] = {{0.0625, 0.125, 0.0625}, {0.125, 0.25, 0.125}, {0.0625, 0.125, 0.0625}};
