#ifndef RAYTRACING_H
#define RAYTRACING_H

#define WINDOW_WIDTH 512
#define WINDOW_HEIGHT 512

#define SAMPLE_RATE 3
extern const double antiAliasMatrix[SAMPLE_RATE][SAMPLE_RATE];

int main(int argc, char *argv[]);

#endif
