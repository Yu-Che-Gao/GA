// GA.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include <time.h>
#include "GA.h"

int main(int argc, char *argv)
{
	int i, j;
	srand((unsigned)time(NULL));
	initialize(); //初始化
	for (i = 0; i < ITERA_CNT; i++)
	{
		reproduction(); //選擇(分配式)
		//reproduction_rnd(); //選擇(隨機式),收斂速度慢
		crossover(); //交配
		mutation(); //突變
		printf("第%d次迭代\n", i);
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

