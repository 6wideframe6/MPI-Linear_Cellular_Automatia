#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <mpi.h>			//���������� ���������� mpi

int TorIt(int x, int count)			//������� �������� ����� ���������� ��������
{
	if (x < 0) return x + count; else return x % count;
}

int f(int y1, int y2, int y3, int choice)				//������� �������� ���������� �������� (�/� XOR)
{
	if (choice == 1)
		return y1 ^ y2 ^ y3;
	if (choice == 2)
		return y1 ^ y3;
	if (choice == 3)
		return y1 | y2 | y3;
}

int numberCount(FILE* input) {				//������� ���������� ����� � �����
	fseek(input, 0, SEEK_SET);
	int counter = 0;
	for (;;) {
		int value;
		if (fscanf(input, "%d", &value) == 1)
			counter++;
		if (feof(input))
			break;
	}
	return counter;
}

void setFirstGen(int* y, int n) {			//������ �� ����� ���������� ������������� ������ ���������
	FILE* fin = NULL;
	fin = fopen("firstgen.txt", "r");
	int i = 0;
	while (i < n)
	{
		fscanf(fin, "%d", &y[i]);
		++i;
	}
	fclose(fin);
	for (i = 0; i < n; i++)
		printf("%d", y[i]);
	printf("\n");
	return y;
}

void shiftMsg(int myid, int numprocs, int n, int steps, int choice) {		//�������� ������� ���������� MPI
	int i;
	int left = n % numprocs;					//������� �� ������� ����� �� ���������� �������� ���������
	int piece = (n - left) / numprocs;			//�������� ��� ������� �������� (����� �������� - ��������)

	int* y = (int*)malloc(sizeof(int) * n);		//��������� ������ � ������ �������� ��� ������ ������ ���������

	int* y1 = NULL;				//���������� ������� ������� ��� ���� ���������
	int length;					//����� ������������ �� ������ ������� �������

	int i1, i2;

	if (myid == 0) {			//�������� �� ������� �������� ������� ������� �������: ������� + �������� ��� ������� ��������
		setFirstGen(y, n);
		y1 = malloc(sizeof(int) * (piece + left));
		length = piece + left;
	}
	else {						//�������� �� ��������� ��������� ������� ������� �������� � �������� ��� ������� ��������
		y1 = malloc(sizeof(int) * (piece));
		length = piece;
		i1 = piece * myid + left;				//��������� ������ ������� ������� ��� ������� ��������, ����� �������� - ��������
		i2 = piece * (myid + 1) + left;
	}


	int* rcounts = NULL;						//��������� �� ������� �������� ������ ���������� ������� ��� MPI_Gatherv
	int* displs = NULL;

	if (myid == 0) {
		rcounts = malloc(numprocs * sizeof(int));
		displs = malloc(numprocs * sizeof(int));
		for (i = 0; i < numprocs; i++) {
			if (i == 0) { 
				rcounts[i] = piece + left; 
				displs[i] = 0;
			}
			else {
				rcounts[i] = piece;
				displs[i] = displs[i-1] + rcounts[i - 1];
			}
		}

	}

	for (int j = 0; j < steps; j++) {				

		MPI_Bcast(y, n, MPI_INT, 0, MPI_COMM_WORLD);	//�������� ����� ������� ���������� ��������� �� ��� ��������

		if (myid == 0)			//������ ������� ������ ��������� ��� �������� ���������� ��������� (y), ��������� ��������� � ���� ������ ������� (y1)
		{
			
			for (int i = 0; i < piece + left; i++)
				y1[i] = f(y[TorIt(i - 1, n)], y[TorIt(i, n)], y[TorIt(i + 1, n)], choice);
		}
		else
		{
			for (int i = i1; i < i2; i++)
				y1[i - i1] = f(y[TorIt(i - 1, n)], y[TorIt(i, n)], y[TorIt(i + 1, n)], choice);
		}

		//���������� ����������� �� ������� ������ � �������� ������ (y) � ������ ������ � �������� ������� ������� (y1) �� ������ ��������
		MPI_Gatherv(y1, length, MPI_INT, y, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

		//����� ���������� �� ������ �������� ����������� ��� ����������� ���������� ������
		if (myid == 0) {
			for (int i = 0; i < n; i++) {
				printf("%d", y[i]);
			}
			printf("\n");
		}

	}

	//������������ ������

	free(y1);

	free(y);

	if (myid == 0) {
		free(rcounts);
		free(displs);
	}

	return;
}

int main(int argc, char* argv[])
{
	int myid, numprocs;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if (argc != 4) {
		if (myid == 0) {
			printf("Needs arg: num of cells, num of iterations and num of method of evolution:\n");
			printf("1: XOR between two neighbour cells and the cell itself\n");
			printf("2: XOR between two neighbour cells (rule 90)\n");
			printf("3: OR between two neighbour cells and the cell itself\n");
		}
		MPI_Finalize();
		return 0;
	}

	int n;						//����������� �����
	n = atoi(argv[1]);			//n - ���������� ������ � �����
	int steps;					//steps - ���������� �������� ��������
	steps = atoi(argv[2]);
	int choice = atoi(argv[3]);

	if (myid == 0) {
		printf("You've chosen:\n");
		printf("Num of cells: %d\n", n);
		printf("Num of iterations: %d\n", steps);
		printf("Method of evolution: %d", choice);
			if (choice == 1)
				printf(" (XOR between two neighbour cells and the cell itself)\n");
			if (choice == 2)
				printf(" (XOR between two neighbour cells (rule 90))\n");
			if (choice == 3)
				printf(" (OR between two neighbour cells and the cell itself)\n");
			printf("\n");
	}

	int error = NULL;

	if (myid == 0) {
		FILE* fin = NULL;
		fin = fopen("firstgen.txt", "r");
		if (!fin) {
			error = 1;
			printf("Error: input file doesnt exist\n");
		}
		if (error == NULL) {
			int num = numberCount(fopen("firstgen.txt", "r"));
			if (n != num) {
				if (n < numberCount(fopen("firstgen.txt", "r"))) {				//�������� ������� ������ � ������� ������ ����� ������������.
					printf("Warning: actual length of input is %d", num);
					printf(", but it will work with num of cells %d", n);
				}
				else {
					printf("Error: actual length of input is %d", num);
					printf(", enter correct num of cells\n");
					error = 1;
				}
			}
			if (choice < 1 || choice > 3) {
				printf("Error: method arg should be from 1 to 3:\n");
				printf("1: XOR between two neighbour cells and the cell itself\n");
				printf("2: XOR between two neighbour cells (rule 90)\n");
				printf("3: OR between two neighbour cells and the cell itself\n");
				error = 1;
			}
		}
		printf("\n");
	}
	MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (error != NULL) {
		MPI_Finalize();
		return 0;
	}

	if (myid == 0) {
		printf("Start: \n");
	}
	shiftMsg(myid, numprocs, n, steps, choice);

	MPI_Finalize();
	return 0;
}
