/*
 ============================================================================
 Name        : LCSSerial.c
 Author      : chathuri
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
				memset(seq1,0,sizeof(seq1));
				strcpy(seq1, line);
			} else if (count == 2) {
				seq2 = (char*) malloc(columnCount);
				memset(seq2,0,sizeof(seq2));
				strcpy(seq2, line);
			} else
				break;
			count++;
		}
	}
	fclose(input);
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

	if(ptContiniousMem==NULL || array==NULL){
		fprintf(stderr,"Can't allocate memory");
	}

	for (i = 0; i <= rowCount; i++)
		array[i] = &(ptContiniousMem[(columnCount + 1) * i]);

	for (i = 1; i <= rowCount; i++) {
		for (j = 1; j <= columnCount; j++) {
			if (seq1[i - 1] == seq2[j - 1]) {
				array[i][j] = array[i - 1][j - 1] + cost(j);
			} else {
				array[i][j] = MAX(array[i - 1][j], array[i][j - 1]);
			}
		}
	}
	/*for(i=0;i<=rowCount;i++){
	 printf("\n");
	 for(j=0;j<=columnCount;j++){
	 fprintf(stderr,"printed %h",array[i][j]);
	 }
	 printf("\n");
	 }*/
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
	memset(commonStr,0,sizeof(commonStr));
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

void clearMemory(){
		free(array[0]);
		free(array);
		free(seq1);
		free(seq2);
}

int main(int argc, char **argv) {

	char *l_pFileName = argv[1];
        double start_t=omp_get_wtime();
	readFromInput(l_pFileName);
	buildCostMatrix();
	traverseToFindMatch();
	double end_t=omp_get_wtime();
	int matchLength=strlen(commonStr);
	fprintf(stdout,"Time taken for serial code: %.16g \n",end_t-start_t);
	fprintf(stdout,"%d\n",matchLength);
	fprintf(stdout, "%s\n", commonStr);
    clearMemory();

	return 0;
}
