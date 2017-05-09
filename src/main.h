#ifndef MAIN_H
#define MAIN_H

#define WINDOW_SIZE 512

#include "input.h"

void display();
void initialize(const std::string &path);
void reshape(int w, int h);
int main(int argc, char *argv[]);

#endif
