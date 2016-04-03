// GA.cpp : �w�q�D���x���ε{�����i�J�I�C
//

#include "stdafx.h"
#include <time.h>
#include "GA.h"

int main(int argc, char *argv)
{
	int i, j;
	srand((unsigned)time(NULL));
	initialize(); //��l��
	for (i = 0; i < ITERA_CNT; i++)
	{
		reproduction(); //���(���t��)
		//reproduction_rnd(); //���(�H����),���ĳt�׺C
		crossover(); //��t
		mutation(); //����
		printf("��%d�����N\n", i);
		for (j = 0; j < POPULATION_CNT; j++)
		{
			printf("(%5.2lf,%5.2lf)", population[j].dec_value, population[j].fitness);
			if (j % 4 == 3) printf("\n");
		}
	}

	printf("\n==================================\n");
	printf("%3d times...\n", i);
	for (j = 0; j < POPULATION_CNT; j++)
	{
		printf("(%5.2lf,%5.2lf)", population[j].dec_value, population[j].fitness);
		if (j % 4 == 3) printf("\n");
	}
	printf("\n==================================\n");
	printf(" ever find best gene:");
	printf("(%5.2lf,%5.2lf)\n", best_gene.dec_value, best_gene.fitness);
	system("pause");
	return 0;
}

