// Minimal Cost of Conversion for a Given Genome
// Author: Henry Ton

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#define g 2.0
#define h 0.5


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

// reverses the sequence with a given length
void reverse(char *sequence, int length)
{
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

float minGapValue(char *sequenceA, char *sequenceB, int N, int aStart)
{
	//if(N == 0)
		//return 0;
	
	float min = gap(0) + w(sequenceA[0], sequenceB[0]) + gap(N);

	for(int j = 1; j <= N; j++)
	{
		if((gap(j - 1) + w(sequenceA[aStart], sequenceB[j]) + gap(N - j)) < min)
			min = (gap(j - 1) + w(sequenceA[aStart], sequenceB[j]) + gap(N - j));
	}
	
	printf("Min gap value: %f\n", min);
	return min;
}

// gap penalty f
double gap(int k)
{
	if (k > 0)
		return g + h*k;
	return 0;
}

// takes in a set of parameters and returns the minimum of those parameters
float min2(float a, float b)
{
	if(a < b) return a;
	
	return b;
}

float min3(float a, float b, float c)
{
	if(a < b && a < c)
		return a;
	else if(b < a && b < c)
		return b;
	else 
		return c;
}

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
	
	printf("Entering alignment calculation\n");
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
	
	// initialize the first index of CC[]
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
			printf("%c - %c\n", sequenceA[i-1], sequenceB[j-1]);
			
			c = min3(array2[j], e, s + w(sequenceA[i-1], sequenceB[j-1]));
			s = array1[j];
			array1[j] = c;
		}
		printf("Cost is %.2f\n", array1[N]);
	}
	array2[0] = array1[0];
	printf("Calculation done\n");
}

float DIFF(char *sequenceA, char *sequenceB, int M, int N)
{
	// diff function that takes a and b indices along with the length of each one M/N
	return diff(sequenceA, sequenceB, 0, 0, M, N, g, g);
}

// recursive procedure for diff
float diff(char *sequenceA, char *sequenceB, int aStart, int bStart, int M, int N, float tb, float te)
{
	
	printf("Calling diff with aStart: %d and bStart: %d\n", aStart, bStart);
	// Denoted as i*, will be the half point index used to spli the array up
	int i, j;
	float first, second, inserted, deleted;
	
	if(N == 0)
	{
		if(M > 0)
		{
			printf("Delete A\n");
			deleted = 1.0;
			return 1.0;
		}
	}
	
	else if(M == 0)
	{
		printf("Insert B\n");
		inserted = 1.0;
		return 1.0;
	}
	
	else if(M == 1)
	{	
		//printf("Conversion of cost min: %f\n", 
			//min(2,
				//(min(2, tb, te) + h) + gap(N), 
					//minGapValue(sequenceA, sequenceB, N, aStart)
		
		float cost1 = (min2(tb, te) + h) + gap(N);
		int type = 1;
		
		for(int j = 1; j <= N; j++)
		{
		
			float cost2 = gap(j - 1) + w(sequenceA[aStart], sequenceB[bStart + j - 1]) + gap(N - 1);
			
			if(cost2 < cost1)
				cost1 = cost2;
		}
		printf("Conversion cost of min: %f\n", cost1);
		return cost1;
	}
	
	else
	{
		// dividing i by 2
		i = M/2;
		printf("dividing i by 2. i* is now %d\n", i);
		
		// Compute CC and DD in a forward phase and reverse phases
		printf("Cumputing CC and DD with M/N being %d/%d\n", i, N);
		alignment('F', i, N, sequenceA, sequenceB, tb);

		reverse(sequenceA, M);
		reverse(sequenceB, N);
		
		printf("Cumputing RR and SS with M/N being %d/%d\n", i, N);
		alignment('R', M - i, N, sequenceA, sequenceB, te);
		
		reverse(sequenceA, M);
		reverse(sequenceB, N);

		// Find minimum value j
		
		j = 0;
		float t1 = CC[0] + RR[N], t2 = DD[0] + SS[N] - g;
		int type = (t1 < t2)? 1: 2;
		float t3 = min2(t1, t2);

		for(int x = 1; x <= N; x++)
		{
			t1 = CC[x] + RR[N - x]; 
			t2 = DD[x] + SS[N - x] - g;
			//t3 = min2(t1, t2);
			float current_min = min2(t1, t2);
			
			if(current_min < t3)
			{
				j = x;
				t3 = current_min;
				type = (t1 < t2)? 1: 2;
			}

		}
		printf("Calculated j value: %d\n", j);
		
		if(type == 1)
		{
			first = diff(sequenceA, sequenceB, aStart, bStart, i, j, tb, g);
			second = diff(sequenceA, sequenceB, aStart + i, bStart + j,  M - i, N - j, g, te);
		}
		else
		{
			first = diff(sequenceA, sequenceB, aStart, bStart, i - 1, j, tb, 0);
			printf("Delete a\n");
			second = diff(sequenceA, sequenceB, aStart + i, bStart + j, M - i - 1, N - j, 0, te);
		}
		
		return first + second;
	}
	
}

void display(char *sequenceA, char  *sequenceB, int M, int N, float *S)
{
	;
}

int main(void) 
{
	// Global variable containing sequences A and B
	// declaring some random DNA sequence
	char sequenceA[5] = {'a', 'a', 'a', 'a', 'a'};
	char sequenceB[4] = {'a', 'g', 'c', 't'};
	
	int M = 5; 
	int N = 4;
	
	char *S = malloc(sizeof(char) * M);
	
	// Allocate memory space to the four arrays
	CC = malloc(sizeof(float) * M + 1);
	DD = malloc(sizeof(float) * M + 1);
	RR = malloc(sizeof(float) * N + 1);
	SS = malloc(sizeof(float) * N + 1);
	S = malloc(sizeof(float) * M + 1);
	
	// Initiate main wrapper for DIFF
	printf("Minimum conversion cost: %f\n" , DIFF(sequenceA, sequenceB, M, N));
	
	free(CC);
	free(DD);
	free(RR);
	free(SS);
	
	return 0;
}

// Matching with a list of people pulled from the database
// edit distance
// deletes are replaced by gaps, 
//  gap{k) = g + hk where g is the fixed cost to open up a gap, h is how much each
// nucleotide costs to remove/delete, k is the number of nucleotides in the gap
// DD is the cost to convert when ending with a deletion
// CC is the cost of conversion - optimal cost of conversion
// I is the cost of conversion when inserting at the last letter
// D is the cost of deletion at the last letter 

// i* is the midpoint of A, j* is the midpoint of B BUT i* doesn't change, j is what you're changing around