/*
 ============================================================================
 Name        : lcs-mpi.c
 Version     :
 Copyright   : Your copyright notice
 Description : Hello MPI World in C 
 ============================================================================
 */
#include <stdio.h>
#include <string.h>
#include<math.h>
#include<stdlib.h>
#include "mpi.h"
int rowCount = 0;
int columnCount = 0;
char *seq1;
char *seq2;
short **array;
short **subArray;
#define MAX(a,b) (((a)>(b))?(a):(b))
char commonStr[30000];

void readFromInput(char *FileName) {
	FILE *input;
	int count = 0;
	char *l_sTem;
	char *tok = NULL;
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
				//printf("row c %d", rowCount);
				tok = strtok(NULL, " ");
				columnCount = atoi(tok);
			} else if (count == 1) {
				seq1 = (char*) malloc(rowCount + 1);
				memset(seq1, 0, sizeof(seq1));
				strcpy(seq1, line);
				seq1[strlen(line) - 1] = '\0';
			} else if (count == 2) {
				seq2 = (char*) malloc(columnCount + 1);
				memset(seq2, 0, sizeof(seq2));
				strcpy(seq2, line);
				seq2[strlen(line) - 1] = '\0';
			} else
				break;
			count++;
		}

	} else {
		printf("no such input");
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

void buildCostMatrix(int displace, int sendCount, int my_rank, int p) {
	short *ptContiniousMem = (short *) calloc((rowCount + 1) * (sendCount),
			sizeof(short));
	subArray = (short **) calloc(rowCount + 1, sizeof(short*));

	if (ptContiniousMem == NULL || subArray == NULL) {
		fprintf(stderr, "Can't allocate memory\n");
	}
	int i, j;

	for (i = 0; i <= rowCount; i++)
		subArray[i] = &(ptContiniousMem[(sendCount) * i]);

	short newBoundary = 0; //last value of same row by previous p
	MPI_Status status;
	MPI_Request *requestList;
	short prevBoundary = 0; //last value of one row above by previous p

	requestList = (MPI_Request*) malloc((rowCount) * sizeof(MPI_Request));
	for (i = 1; i <= rowCount; i++) {
		if (displace != 0) {
			MPI_Recv(&newBoundary, 1, MPI_SHORT, my_rank - 1, i, MPI_COMM_WORLD,
					&status);
		}
		for (j = 0; j < sendCount; j++) {
			if (displace == 0 && j == 0) {
				continue;
			}

			if (seq1[i - 1] == seq2[displace + j - 1]) {
				if (j == 0) {
					subArray[i][j] = prevBoundary + cost(j);
				} else {
					subArray[i][j] = subArray[i - 1][j - 1] + cost(j);
				}
			} else {
				if (j == 0) {
					subArray[i][j] = MAX(subArray[i - 1][j], newBoundary);
				} else {
					subArray[i][j] = MAX(subArray[i - 1][j],
							subArray[i][j - 1]);
				}
			}
		}

		if (displace != 0) {
			prevBoundary = newBoundary;
		}
		newBoundary = subArray[i][sendCount - 1];
		if (my_rank < p - 1) {
			MPI_Isend(&newBoundary, 1, MPI_SHORT, my_rank + 1, i,
					MPI_COMM_WORLD, &(requestList[i - 1]));
		}
	}

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

//each processor will find their local match and at the end root will do a gatherv
void traverseToFindMatch(int numColumns, int displace, int my_rank, int p,
		int commonStrCnt[]) {
	memset(commonStr, 0, sizeof(commonStr));
	MPI_Status status;
	MPI_Status statusArr;
	MPI_Request request;
	int localEnd = numColumns - 1;
	short currentCell;
	int i = rowCount;
	int j = localEnd;
	int count = 0;
	int newPos;
	int globalColumnPos = displace + numColumns - 1;

	memset(commonStr, 0, sizeof(commonStr));
	short *boundaryArray = (short*) malloc((rowCount + 1) * sizeof(short));
	short *receiveArray = (short*) malloc((rowCount + 1) * sizeof(short));

	if (my_rank < p - 1) {
		//send the boundary
		int k;
		for (k = 0; k <= rowCount; k++) {
			boundaryArray[k] = subArray[k][localEnd];
		}

		MPI_Send(boundaryArray, rowCount + 1, MPI_SHORT, my_rank + 1, 000,
				MPI_COMM_WORLD);
		MPI_Recv(&newPos, 1, MPI_INT, my_rank + 1, 010, MPI_COMM_WORLD,
				&status);

		i = newPos;
		currentCell = (i == 0) ? 0 : subArray[i][localEnd];
	} else {
		currentCell = subArray[rowCount][localEnd];
	}
	if (my_rank > 0) {
		MPI_Irecv(receiveArray, rowCount + 1, MPI_SHORT, my_rank - 1, 000,
				MPI_COMM_WORLD, &request); //non blocking as we need this array only at the end
	}
	while (currentCell > 0 && i > 0) {
		if (j > 0) {
			if (seq1[i - 1] == seq2[globalColumnPos - 1]) {
				i = i - 1;
				j = j - 1;
				globalColumnPos = globalColumnPos - 1;
				char let = seq1[i];
				commonStr[count] = let;
				count++;
			} else {
				if (subArray[i][j - 1] >= subArray[i - 1][j]) {
					j = j - 1;
					globalColumnPos = globalColumnPos - 1;
				} else {
					i = i - 1;
				}
			}
			currentCell = subArray[i][j];
		}

		else {
			MPI_Wait(&request, &statusArr);	//u need the receive buffer after this.
			if (currentCell == 0 && my_rank > 0) {
				newPos = 0;
				MPI_Send(&newPos, 1, MPI_INT, my_rank - 1, 010, MPI_COMM_WORLD);
				break;
			}
			if (seq1[i - 1] == seq2[globalColumnPos - 1]) {
				i = i - 1;
				newPos = i;
				MPI_Send(&newPos, 1, MPI_INT, my_rank - 1, 010, MPI_COMM_WORLD);
				char let = seq1[i];
				commonStr[count] = let;
				count++;
				break;
			} else {
				if (receiveArray[i] >= subArray[i - 1][j]) {
					newPos = i;
					MPI_Send(&newPos, 1, MPI_INT, my_rank - 1, 010,
							MPI_COMM_WORLD);
					break;
				} else {
					i = i - 1;
					currentCell = subArray[i][j];
				}
			}
		}
	}
	commonStr[count] = '\0';
	int matchLength = strlen(commonStr);
	commonStrCnt[my_rank] = matchLength;
	int *ptCountRemoteArray;
	ptCountRemoteArray = (int*) malloc(1 * sizeof(int));
	ptCountRemoteArray[0] = commonStrCnt[my_rank];
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Gather((void*) ptCountRemoteArray, 1, MPI_INT, (void*) commonStrCnt, 1,
			MPI_INT, 0, MPI_COMM_WORLD);
	Reverse(commonStr);
}

void clearMemory() {
	free(subArray[0]);
	free(subArray);
	free(seq1);
	free(seq2);
}

int main(int argc, char* argv[]) {
	int my_rank; /* rank of process */
	int p; /* number of processes */
	int source; /* rank of sender */
	MPI_Status status; /* return status for receive */
	double elapsedTime;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &p);

	elapsedTime -= MPI_Wtime();
	int i;
	if (my_rank == 0) { //main process is going to load the data
		char *l_pFileName = argv[1];
		readFromInput(l_pFileName);
	}

	MPI_Bcast(&rowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&columnCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (my_rank != 0) {
		seq1 = (char*) malloc(rowCount + 1);
		memset(seq1, 0, sizeof(seq1));
		seq2 = (char*) malloc(columnCount + 1);
		memset(seq2, 0, sizeof(seq2));
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(seq1, rowCount + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(seq2, columnCount + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

	int firstElement = (my_rank * (columnCount + 1)) / p;
	int lastElement = (((my_rank + 1) * (columnCount + 1)) / p) - 1;
	int sendCount = lastElement + 1 - firstElement;
	buildCostMatrix(firstElement, sendCount, my_rank, p);
	MPI_Barrier(MPI_COMM_WORLD);

	//the size of matching string each process is calculated and will be gathered
	int commonStrCnt[p];
	int commonStrDisplacement[p];
	traverseToFindMatch(sendCount, firstElement, my_rank, p, commonStrCnt);

	char *combinedCommonStr;
	if (my_rank == 0) {
		int k = 0;
		for (i = 0; i < p; i++) {
			commonStrDisplacement[i] = k;
			k = k + commonStrCnt[i];
		}
		combinedCommonStr = (char*) malloc(k + 1);
		memset(combinedCommonStr, 0, sizeof(k + 1));
	}
	//Process 1 is going to do a gather to create combined string
	MPI_Gatherv((void*) commonStr, commonStrCnt[my_rank], MPI_CHAR,
			(void*) combinedCommonStr, commonStrCnt, commonStrDisplacement,
			MPI_CHAR, 0, MPI_COMM_WORLD);

	elapsedTime += MPI_Wtime();
	if (my_rank == 0) {
		combinedCommonStr[strlen(combinedCommonStr)] = '\0';
		int length = strlen(combinedCommonStr);
		printf("%d\n", length);
		printf("%s\n", combinedCommonStr);
		fprintf(stderr,"time difference is %10.6f\n", elapsedTime);
	}
	clearMemory();
	MPI_Finalize();

	return 0;
}
