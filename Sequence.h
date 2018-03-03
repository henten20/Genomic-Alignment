#ifndef SEQUENCE
#define SEQUENCE

#define g 2.0
#define h 0.5
#define MAX_MEM_USAGE 1000

// globally shared vectors for forward and reverse phases
float *CC;
float *DD;
float *RR;
float *SS;
char *s;

// Functional Prototypes
double gap(int);
float w(char, char);
void alignment(char, int, int, char *, char *, float);
float DIFF(char *, char *, int, int);
float diff(char *, char *, int, int, int, int, float, float);

#endif