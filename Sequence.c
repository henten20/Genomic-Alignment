// Minimal Cost of Conversion for a Given Genome
// Author: Henry Ton

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "Sequence.h"

// reverses the sequence with a given length
void reverse(char *sequence, int length)
{
	if(sequence == NULL || length < 1)
		return;
		
	char temp;
	int i = 0, j = length - 1;
	
	while(i < j)
	{
		temp = sequence[i];
		sequence[i] = sequence[j];
		sequence[j] = temp;
		
		i++; j--;
	}
}

// gap penalty f
double gap(int k)
{
	if (k > 0)
		return g + h*k;
	
	return 0;
}

// return minimum of two given arguments
float min2(float a, float b)
{
	if(a < b) return a;
	
	return b;
}

// return minimum of three given arguments
float min3(float a, float b, float c)
{
	if(a < b && a < c)
		return a;
	else if(b < a && b < c)
		return b;
	else 
		return c;
}

// return the cost to convert a to b
// 0 if they are equal.
float w(char a, char b)
{
	if(a == b)
		return 0.0;
	
	else return 1.0;
}

// calculate the minimal cost of conversion for a sub sequence
void alignment(char phase, int M, int N, char *sequenceA, char *sequenceB, float x)
{
	float t;
	float *array1, *array2;
	
	// 'F' or 'R' dictates a forward or reverse phase
	if(phase == 'F')
	{
		array1 = CC;
		array2 = DD;
	}
	else
	{
		array1 = RR;
		array2 = SS;
	}
	
	array1[0] = 0;
	t = g;
	
	for(int j = 1; j <= N; j++)
	{
		t += h;
		array1[j] = t;
		array2[j] = t + g;
	}
	
	// [*]
	t = x;

	for(int i = 1; i <= M; i++)
	{
		float s, c, e;
		s = array1[0];
		t += h;
		c = t;
		array1[0] = c;
		e = t + g;
	
		for(int j = 1; j <= N; j++)
		{
			e = min2(e, c + g) + h;
			array2[j] = min2(array2[j], array1[j] + g) + h;
			c = min3(array2[j], e, s + w(sequenceA[i-1], sequenceB[j-1]));
			s = array1[j];
			array1[j] = c;
		}
		
		//printf("Cost is %.2f\n", array1[N]);
	}
	
	array2[0] = array1[0];
}

float DIFF(char *sequenceA, char *sequenceB, int M, int N)
{
	// diff function that takes a and b indices along with the length of each one M/N
	return diff(sequenceA, sequenceB, 0, 0, M, N, g, g);
}

// recursive procedure for diff
float diff(char *sequenceA, char *sequenceB, int aStart, int bStart, int M, int N, float tb, float te)
{
	if(sequenceA == NULL || sequenceB == NULL)
		return -1.0;
	
	int i, j;
	float forwardAlignment, reverseAlignment;
	
	if(N == 0)
	{
		if(M > 0)
		{
			//printf("Delete A\n");
			return 1.0;
		}
	}
	
	else if(M == 0)
	{
		//printf("Insert B\n");
		return 1.0;
	}
	
	else if(M == 1)
	{	
		float cost1 = (min2(tb, te) + h) + gap(N);
		int type = 1;
		
		for(int j = 1; j <= N; j++)
		{
			float cost2 = gap(j - 1) + w(sequenceA[aStart], sequenceB[bStart + j - 1]) + gap(N - 1);
			
			if(cost2 < cost1)
				cost1 = cost2;
		}
		
		//printf("Conversion of cost min is %f\n", cost1);
		return cost1;
	}
	
	else
	{
		i = M/2;
		
		// Compute CC and DD in a forward phase
		alignment('F', i, N, sequenceA, sequenceB, tb);

		reverse(sequenceA, M);
		reverse(sequenceB, N);
		
		// Compute RR and SS in a reverse phase
		alignment('R', M - i, N, sequenceA, sequenceB, te);
		
		reverse(sequenceA, M);
		reverse(sequenceB, N);

		// Find minimum value j
		j = 0;
		
		float t1 = CC[0] + RR[N];
		float t2 = DD[0] + SS[N] - g;
		float t3 = min2(t1, t2);
		int type = (t1 < t2)? 1: 2;
		
		for(int x = 1; x <= N; x++)
		{
			t1 = CC[x] + RR[N - x]; 
			t2 = DD[x] + SS[N - x] - g;

			float current_min = min2(t1, t2);
			
			if(current_min < t3)
			{
				j = x;
				t3 = current_min;
				type = (t1 < t2)? 1: 2;
			}
		}
		
		if(type == 1)
		{
			forwardAlignment = diff(sequenceA, sequenceB, aStart, bStart, i, j, tb, g);
			reverseAlignment = diff(sequenceA, sequenceB, aStart + i, bStart + j,  M - i, N - j, g, te);
		}
		else
		{
			forwardAlignment = diff(sequenceA, sequenceB, aStart, bStart, i - 1, j, tb, 0);
			//printf("Delete A\n");
			reverseAlignment = diff(sequenceA, sequenceB, aStart + i, bStart + j, M - i - 1, N - j, 0, te);
		}
		
		return forwardAlignment + reverseAlignment;
	}
}

int main(int argc, char *argv[]) 
{
	if(argc < 3)
	{
		printf("Incorrect format.\n");
		return 0;
	}

	FILE *A, *B;
	A = fopen(argv[1], "r");
	B = fopen(argv[2], "r");
	
	if(A == NULL || B == NULL)
	{
		printf("Error in locating sequences.\n");
		exit(1);
	}

	int M = 0, N = 0;
	char seqA[MAX_MEM_USAGE], seqB[MAX_MEM_USAGE];
	
	fgets(seqA, MAX_MEM_USAGE, A);
	fgets(seqB, MAX_MEM_USAGE, B);
	
	printf("%s\n", seqA);
	printf("%s\n", seqB);
	
	for(int i = 0; ;i++)
	{
		if(seqA[i] == '\0')
			break;
		M++;
	}
	
	for(int i = 0; ;i++)
	{
		if(seqB[i] == '\0')
			break;
		N++;
	}
	
	char *sequenceA = malloc(sizeof(char) * M);
	char *sequenceB = malloc(sizeof(char) * N);
	
	for(int i = 0; i < M; i++)
		sequenceA[i] = seqA[i];
	
	for(int i = 0; i < N; i++)
		sequenceB[i] = seqB[i];
	
	//char sequenceA[7] = {'a', 'a', 'a', 'a', 'a', 'a', 't'};
	//char sequenceB[5] = {'a', 'a', 'c', 't', 'c'};
	//char *S = malloc(sizeof(char) * M);
	
	// Allocate memory space to the four arrays
	CC = malloc(sizeof(float) * M + 1);
	DD = malloc(sizeof(float) * M + 1);
	RR = malloc(sizeof(float) * N + 1);
	SS = malloc(sizeof(float) * N + 1);
	//S = malloc(sizeof(float) * M + 1);
	
	printf("Minimum conversion cost: %.2f\n" , DIFF(sequenceA, sequenceB, M, N));
	
	free(sequenceA);
	free(sequenceB);
	free(CC);
	free(DD);
	free(RR);
	free(SS);
	
	fclose(A);
	fclose(B);
	
	return 0;
}
