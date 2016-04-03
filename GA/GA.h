#pragma once
#ifndef __GA__
#define __GA__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define GENETIC_LENGTH 4 //��]����
#define POPULATION_CNT 10 //���s�ƶq
#define ITERA_CNT 100 //���N����
#define CROSSOVER_RATE 0.2 //��t�v
#define MUTATION_RATE 0.1 //���ܲv

//--------------------------------------
//�w�q���鵲�c
typedef struct parent_t {
	int genes[GENETIC_LENGTH];
	double fitness;
	double dec_value;
}praent_t;
//--------------------------------------

//GAPosRand(): �H�����o���ܦ�m
#define GAPosRand() (rand()%GENETIC_LENGTH)

//BinaryRand(): �H������/1 �����
#define BinaryRand() (rand()%2)

//SRand(): �H������~1�����
#define SRand() ((double)rand()/(double)RAND_MAX)


//-----------------------------------------
void initialize(); //�i���l��
void reproduction(); //�ƻs,���L�����(���t��),�M�w�C�ӥ���ƻs���t�����Ӽ�
void reproduction_rnd(); //�ƻs,���L�����(�H����)
void crossover(); //��t,��t�����������t,[���I��t,���I��t,mask]
void mutation(); //����,�v�@bit�C�C�T�{����
void cal_fitness(parent_t *x); //�p���]�ҹ������A����
void cal_xvalue(parent_t *x); //�p���]�������i���
//-----------------------------------------

parent_t population[POPULATION_CNT]; //����ƶq
parent_t pool[POPULATION_CNT]; //��t��
parent_t best_gene; //�q�H�e��{�b�̦n����]

//------------------------------------------
//binary 2 dec,�N�V���餤��genes�ഫ���Q�i��
void cal_xvalue(parent_t *x)
{
	int i, dec = 0;
	for (i = 0; i < GENETIC_LENGTH; i++)
	{
		if (x->genes[i] == 1) dec = dec + (0x01 << i);
	}
	x->dec_value = (double)dec;
}
//------------------------------------------

//------------------------------------------
//�A���禡,���]��f(x)=x*x, x���V���餤���Q�i��, �Ydec_value
void cal_fitness(parent_t *x)
{
	double i = x->dec_value;
	x->fitness = (-1)*i*i;
}
//------------------------------------------

//------------------------------------------
//��l��
void initialize()
{
	int i, j;
	for (i = 0; i < POPULATION_CNT; i++)
	{
		for (j = 0; j < GENETIC_LENGTH; j++)
		{
			population[i].genes[j] = BinaryRand(); //�C�ӥ��鳣�O�H���� /1
		}
		cal_xvalue(&population[i]); //�p������]���i���
		cal_fitness(&population[i]); //�p�����������A����

		if (i == 0)
		{
			memcpy(&best_gene, &population[i], sizeof(parent_t));
		}
		else if (population[i].fitness>best_gene.fitness)
		{
			memcpy(&best_gene, &population[i], sizeof(parent_t));
		}
	}
}
//-------------------------------------------

//-------------------------------------------
//�ƻs,���L�����(���t��)
void reproduction()
{
	int i, j, cnt, has_copy = 0;
	int Slack = POPULATION_CNT;//�ٳѴX�ӥi�ƻs
	int pos1, pos2;
	double fitness_sum = 0.0;
	for (int i = 0; i < POPULATION_CNT; i++) //�p��Ҧ��A�����`�M
	{
		fitness_sum += population[i].fitness;
	}

	for (i = 0; i < POPULATION_CNT && Slack != 0; i++) //�p��C�ӥ������ƻs�X�Ө��t����,�ê����@�ƻs
	{
		cnt = (int)(population[i].fitness / fitness_sum + 0.5); //�p��ƻs�Ӽ�,�|�ˤ��J
		if (cnt > Slack) cnt = Slack;
		for (j = 0; j < cnt; ++j, ++has_copy)
		{
			memcpy(&best_gene, &population[i], sizeof(parent_t));
		}
		Slack -= cnt;
	}

	while (has_copy < POPULATION_CNT) //�Y�٦��S�ƻs����
	{
		pos1 = rand() % POPULATION_CNT; //�H���D������P���V����X��
		do
		{
			pos2 = rand() % POPULATION_CNT;
		} while (pos1 == pos2);
		if (population[pos1].fitness>population[pos2].fitness) i = pos1; //����n����������t��
		memcpy(&pool[has_copy++], &population[i], sizeof(parent_t));
	}
}
//--------------------------------------------

//--------------------------------------------
//�ƻs,���L�����(�H����)
void reproduction_rnd()
{
	int i, pos;
	double fitness_sum = 0.0; //�A�����`�M
	double column_prob[POPULATION_CNT]; //�֭p���v
	double prob; //�H�����v

	for (i = 0; i < POPULATION_CNT; i++)
	{
		fitness_sum += population[i].fitness;
	}
	column_prob[0] = population[0].fitness / fitness_sum;
	for (i = 0; i < POPULATION_CNT; ++i)
	{
		column_prob[i] = column_prob[i - 1] + population[i].fitness / fitness_sum;
	}

	for (i = 0; i < POPULATION_CNT; ++i)
	{
		prob = SRand(); //���Ͷü�
		for (pos = 0; pos < POPULATION_CNT; ++pos)
		{
			if (prob >= column_prob[pos]) break;
		}
		memcpy(&pool[i], &population[pos], sizeof(parent_t));
	}
}
//--------------------------------------------

//--------------------------------------------
//��t
void crossover()
{
	int i, itera;
	int cnt = 0;
	int pos = 0;
	int p1, p2;
	double crossover_if;

	for (itera = 0; itera < POPULATION_CNT; itera++)
	{
		p1 = rand() % POPULATION_CNT;//�H�����ӭ���
		do
		{
			p2 = rand() % POPULATION_CNT;
		} while (p2 == p1);

		crossover_if = SRand(); //�M�w�O�_��t
		if (crossover_if > CROSSOVER_RATE)
		{
			memcpy((void *)&population[cnt++], (void *)&pool[p1], sizeof(parent_t));
			memcpy((void *)&population[cnt++], (void *)&pool[p2], sizeof(parent_t));
		}
		else
		{
			do
			{
				pos = GAPosRand(); //���I��t,��t�����^����
			} while (pos == 0);

			for (i = 0; i < pos; i++) //crossover
			{
				population[cnt].genes[i] = pool[p1].genes[i];
				population[cnt + 1].genes[i] = pool[p2].genes[i];
			}
			cnt += 2; //�w�ƻs�����
		}

	}
}
//--------------------------------------------

//--------------------------------------------
//����
void mutation()
{
	int i;
	int pos;
	for (i = 0; i < POPULATION_CNT; i++)
	{
		double mutation_if = SRand();
		if (mutation_if <= MUTATION_RATE)
		{
			pos = GAPosRand(); //���ܦ�m
			population[i].genes[pos] = 1 - population[i].genes[pos];
		}

		//���ܧ���A��@������A����
		cal_xvalue(&population[i]); //���p���]������x��
		cal_fitness(&population[i]); //�A�N�i��x�ȱa�J�A���禡
		//�A��sbest_gene
		if (i == 0)
		{
			memcpy(&best_gene, &population[i], sizeof(parent_t));
		}
		else if (population[i].fitness>best_gene.fitness)
		{
			memcpy(&best_gene, &population[i], sizeof(parent_t));
		}
	}
}
//--------------------------------------------

#endif // !__GA__