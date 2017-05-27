#ifndef UPDATION_H
#define UPDATION_H

const int MAXNLINE = 1000;
const int MAXFORKS = 4;

void sig_pipe(int signo);
void listenTerminal();
void loadHitpoints();
void updateHitpoints();
int main(int argc, char* argv[]);

#endif
