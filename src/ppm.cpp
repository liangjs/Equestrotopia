#include "ppm.h"

#if (SAMPLE_RATE == 3)
const double antiAliasMatrix[SAMPLE_RATE][SAMPLE_RATE] = {{0.0625, 0.125, 0.0625}, {0.125, 0.25, 0.125}, {0.0625, 0.125, 0.0625}};
#endif
#if (SAMPLE_RATE == 2)
const double antiAliasMatrix[SAMPLE_RATE][SAMPLE_RATE] = {{0.25, 0.25}, {0.25, 0.25}};
#endif
