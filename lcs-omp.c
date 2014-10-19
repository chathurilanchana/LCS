/*
 ============================================================================
 Name        : lcs-omp.c
 Author      : Ashansa, Shelan, Chathuri
 Version     :
 Copyright   : Your copyright notice
 Description : Hello OpenMP World in C
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define MAX(a,b) (((a)>(b))?(a):(b))

int rowCount;
int columnCount;
char *seq1;
char *seq2;
short **array;
char commonStr[30000];

void readFromInput(char *FileName) {
	FILE *input;
	int count = 0;
	char *l_sTem;
	char *tok = NULL;
	//printf(_pFileName);
	input = fopen(FileName, "r");
	if (input != NULL) {

		char line[50000];
		while (fgets(line, sizeof line, input) != NULL) /* read a line */
		{

			if (count == 0) {
				l_sTem = (char*) malloc(100);
				strcpy(l_sTem, line);
				tok = strtok(l_sTem, " ");
				rowCount = atoi(tok);
				tok = strtok(NULL, " ");
				columnCount = atoi(tok);
			} else if (count == 1) {
				seq1 = (char*) malloc(rowCount);
				memset(seq1, 0, sizeof(seq1));
				strcpy(seq1, line);
			} else if (count == 2) {
				seq2 = (char*) malloc(columnCount);
				memset(seq2, 0, sizeof(seq2));
				strcpy(seq2, line);
			} else
				break;
			count++;
		}
	}
	fclose(input);
	free(input);
}

short cost(int x) {
	int i, n_iter = 20;
	double dcost = 0;
	for (i = 0; i < n_iter; i++)
		dcost += pow(sin((double) x), 2) + pow(cos((double) x), 2);
	return (short) (dcost / n_iter + 0.1);

}

void buildCostMatrix() {
	int i, j;

	short *ptContiniousMem = (short *) calloc(
			(rowCount + 1) * (columnCount + 1), sizeof(short));
	array = (short **) calloc(rowCount + 1, sizeof(short*));

	if (ptContiniousMem == NULL || array == NULL) {
		fprintf(stderr, "Can't allocate memory");
	}

  #pragma omp parallel for private(i)
	for (i = 0; i <= rowCount; i++)
		array[i] = &(ptContiniousMem[(columnCount + 1) * i]);

	int hasLengthMaxReached = 0;
	int innerLoopStartingIndex = 1;
	int innerLoopLimit = 0;
	int previousInnerLoopLimit = 0;
	int overflowcounter = 1;
	for (i = 1; i <= rowCount;) {

		innerLoopLimit = i;
		if (hasLengthMaxReached == 1) {
			innerLoopStartingIndex = overflowcounter;
			if (rowCount < columnCount) {
				innerLoopLimit =
						(previousInnerLoopLimit + 1) > columnCount ?
								columnCount : (previousInnerLoopLimit + 1);
			}
		}
		if (i > columnCount) {
			innerLoopLimit = columnCount;
		}

		int y;

		/*We have removed the dependency in the loop by calculating values
		 diagonally. We have created a private copy of val y for each thread and initialize it
		 inside the for loop for each thread. If we use x-- for the second line we have to use locks
		 and have to avoid races. But with this,we dont need them*/

	//the default scheduling which is chunk size=(loop count)/num_thread is the best for us as each iteraton has the same workload
   #pragma omp parallel for private(j,y) if((innerLoopLimit-innerLoopStartingIndex)>=100) //schedule (static,100)
		for (j = innerLoopStartingIndex; j <= innerLoopLimit; j++) {
			y = i - (j - overflowcounter);
			if (seq1[y - 1] == seq2[j - 1]) {
				array[y][j] = array[y - 1][j - 1] + cost(j);
			} else {
				array[y][j] = MAX(array[y - 1][j], array[y][j - 1]);
			}
		}

		previousInnerLoopLimit = innerLoopLimit;
		if (i == rowCount) {
			hasLengthMaxReached = 1;
			overflowcounter++;
			if (overflowcounter > columnCount)
				break;
		} else {
			i++;
		}

	}
	/*for(i=0;i<=rowCount;i++){
	 printf("\n");
	 for(j=0;j<=columnCount;j++){
	 fprintf(stderr,"printed %h",array[i][j]);
	 }
	 printf("\n");
	 }  */
}

void Reverse(char str[]) {
	int i, j;
	char temp[30000];
	for (i = strlen(str) - 1, j = 0; i + 1 != 0; --i, ++j) {
		temp[j] = str[i];
	}
	temp[j] = '\0';
	strcpy(str, temp);
}

void traverseToFindMatch() {

	short currentCell = array[rowCount][columnCount];
	int i = rowCount;
	int j = columnCount;
	memset(commonStr, 0, sizeof(commonStr));
	int count = 0;
	while (currentCell > 0) {
		if (seq1[i - 1] == seq2[j - 1]) {
			i = i - 1;
			j = j - 1;
			char let = seq1[i];
			commonStr[count] = let;

			count++;
		} else {
			if (array[i][j - 1] > array[i - 1][j]) {
				j = j - 1;
			} else if (array[i][j - 1] < array[i - 1][j]) {
				i = i - 1;
			} else {
				j = j - 1;
			}
		}

		currentCell = array[i][j];

	}
	commonStr[count] = '\0';
	//fprintf(stderr, "string is %s", commonStr);
	Reverse(commonStr);
}

void clearMemory() {
	free(array[0]);
	free(array);
	free(seq1);
	free(seq2);
}

int main(int argc, char **argv) {

	/* used to test the correctness of program by setting this to more than cores
#pragma omp parallel
	{
	omp_set_num_threads(10);
    }*/

	char *l_pFileName = argv[1];
	readFromInput(l_pFileName);
	time_t start_t, end_t;
	time(&start_t);
	buildCostMatrix();
	traverseToFindMatch();
	time(&end_t);
	double diff_t;
	diff_t = difftime(end_t, start_t);
	int matchLength = strlen(commonStr);
	/*fprintf(stdout, "Time taken for Building the matrix and find LCS=%f\n",
			diff_t);*/
	fprintf(stdout, "%d\n", matchLength);
	fprintf(stdout, "%s\n", commonStr);
	clearMemory();

	return 0;
}

