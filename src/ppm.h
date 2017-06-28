#ifndef PPM_H
#define PPM_H

#define WINDOW_WIDTH 512
#define WINDOW_HEIGHT 512
#define SAMPLE_RATE 2
extern const double antiAliasMatrix[SAMPLE_RATE][SAMPLE_RATE];

#define PHOTONSPER 500000

#define MAXBOUNCES 30

#define MAXITERATION 15

#define INITRADIUS_RATIO 0.45

#define ALPHA 0.65

#define RUNTHREADS 7

#define PHOTON_CONTRIB 20

#endif
