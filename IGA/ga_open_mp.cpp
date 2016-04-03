#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>
#include<algorithm>
#include<time.h>
#include <string>
using namespace std;
#include <iostream>
#include<sstream>
#include <fstream>

#define M_PI 3.14159265358979323846
#define RANGE_LOWER 0
#define RANGE_UPPER 1

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

static const int population_size=30;//�s���
static const int number_of_gene=4;//�����


static const int max_generation=5000;//�|�N��
float cross_over_probability=0.6f;//��t�v
float alpha=0.25f;//alpha��
float mutation_probability=0.1f;//���ܲv
bool elitism=true;//��׭^�O�d
float scaling_constant=1.2f;//�u���ഫ
int count_times=1;//�D�Ѧ���
int terminate_generation=500;//����@�N��
int restgenerate=0;//�p�⤤��@�N����
//�h�׭^�O�d
float keepfitness[3];
float keepunknowns[3][number_of_gene];
int goodgeneration=0;

#include"store.cuh"
FILE *fp;
#include"osa.cuh"


float generate_rand_zero_to_one()//����0~1��������
{
	return (float)rand()*(1.0f / RAND_MAX);
}

void reset()//�ܼƪ���l��
{
#pragma omp parallel for
	for (individual = 0; individual < population_size; individual++)
	{
		for (gene = 0; gene < number_of_gene; gene++)
		{
			unknowns[individual][gene] = 0.0f;
		}
	}


#pragma omp parallel for
	for (individual = 0; individual < population_size; individual++)
	{
		fitness[individual] = 0.0f;
	}

#pragma omp parallel for
	for (individual = 0; individual < population_size; individual++)
	{
		for (gene = 0; gene < number_of_gene; gene++)
		{
			new_unknowns[individual][gene] = 0.0f;
		}
	}

#pragma omp parallel for
	for (gene = 0; gene < number_of_gene; gene++)
	{
		elite_unknowns[gene] = 0.0f;
	}

	elite_fitness = -1000000000.0f;

	mean_fitness = 0;
	sum_fitness = 0;
	fittest_individual = 0;
	mate1 = 0;
	mate2 = 0;
	generation = 0;
	new_individual = 0;

}

void define_range(float lower, float upper)//�w�q��]���W�U��
{
//#pragma omp parallel for private(gene)
//	for (gene = 0; gene < number_of_gene; gene++)
//	{
//		range[RANGE_LOWER][gene] = lower;
//		range[RANGE_UPPER][gene] = upper;
//	}
}

void initial_population()//��l�Ƹs��
{
#pragma omp parallel
	{
		srand((unsigned int)time(NULL)^omp_get_thread_num());
	#pragma omp for
		for (individual = 0; individual < population_size; individual++)//�H����X�h�ӹ��
		{
			for (gene = 0; gene < number_of_gene; gene++)
			{
				unknowns[individual][gene] = RANGE_LOWER + generate_rand_zero_to_one()*(RANGE_UPPER - RANGE_LOWER);
			}
		}
	}
}


void calculate_fitness()//�p��A����
{
		xsum = 0, ysum = 0 , xysum = 0, x2sum = 0, y2sum = 0;
		for (individual = 0; individual < population_size; individual++)
		{
				fitness[individual] = 0;

				BL_time = (int)(90 * unknowns[individual][0] + 31);
				BL_pre = (int)(80 * unknowns[individual][1] + 6);
				Decline = (int)(3 * unknowns[individual][2] + 2);
				continued = (int)(5 * unknowns[individual][3] + 2);
				
#pragma omp parallel
	{
		BLwm = (float*)malloc(sizeof(float) * BL_time);
		BLarr = (float*)malloc(sizeof(float) * (ROW-BL_time));
	#pragma omp for
				for(int p=0;p<datas;p++) {		
						for(int i=0; i!=BL_time; ++i) {
							BLwm[i] = 0;
						}
						for(int i=0; i!=(ROW-BL_time); ++i) {
							BLarr[i] = 0;
						}
							times=0,ODI=0,Judge=1;


					BL_top = (int)ceil(BL_time*((float)BL_pre/100));

					for(int e = BL_time ; e < ROW2[p] ; e++) {

						for(int m=e-BL_time ; m < e ; m++) {
							BLwm[e - m] = arri[p][m];
						}     
						sort(BLwm, BLwm+BL_time);
				
						top_mean = 0;
						for(int k=0 ; k < BL_top ; k++) { 
							top_mean += (BLwm[k + ((BL_time - 1) - BL_top)]);         
						}
						BLarr[e-BL_time] = top_mean / (float)BL_top;
					}


					for(int i=BL_time;i<ROW2[p];i++) {
						if((((arri[p][i]/BLarr[i-BL_time])-1.0000)*(-100))>=Decline)
						{
							times++;
						}
						else{ 
							times=0;
							Judge=1;
						}
						if((times>=continued))
						{
							if(Judge==1)
							{
								ODI++;
								Judge=0;
							}
						}
					}

					if(ROW2[p] != 0) {
						ODI_L[p]=(float)ODI/((float)ROW2[p]/3600);
					} else {
						ODI_L[p]=0;
					}
				}
	}

					xsum = 0, ysum = 0 , xysum = 0, x2sum = 0, y2sum = 0;
				
					for(int i=0;i< datas;i++) {
						xsum += ODI_L[i];
					}
					for(int i=0;i< datas;i++) {
						ysum += AHI[i];
					}

					xsum = xsum / ((float)datas);
					ysum = ysum / ((float)datas);

					for(int i=1;i< (int)datas;i++) {
							xysum += (ODI_L[i] - xsum) * (AHI[i] - ysum);
							x2sum += pow((ODI_L[i] - xsum), 2);
							y2sum += pow((AHI[i] - ysum), 2);
					}
						fitness[individual] = (xysum / (sqrt(x2sum) * sqrt(y2sum)));
		}
}


void statistics()//�׭^���N
{
	float max_fitness = -1000000000.0f;
	for (individual = 0; individual < population_size; individual++)
	{
		if (fitness[individual]>max_fitness)
		{
			max_fitness = fitness[individual];
			fittest_individual = individual;
		}
	}

	if (elitism == true)//��׭^�O�d
	{
		if (fitness[fittest_individual] <= elite_fitness)//�p�G�ϰ�̨θѥ��j�����̨θ�
		{
				individual = (int)floorf((float)(population_size - 1)*(float)(generate_rand_zero_to_one()));
				fitness[individual] = elite_fitness;
#pragma omp parallel for
				for (gene = 0; gene < number_of_gene; gene++)
				{
					unknowns[individual][gene] = elite_unknowns[gene];
				}
				fittest_individual = individual;
				
//�����|�N
			restgenerate += 1;
			goodgeneration += 1;
		}
		else
		{
			restgenerate = 0;
			goodgeneration = 0;
		}

		elite_fitness = fitness[fittest_individual];//�O�s�׭^�A����ƭ�

#pragma omp parallel for
		for (gene = 0; gene < number_of_gene; gene++)
		{
			elite_unknowns[gene] = unknowns[fittest_individual][gene];//�O�s�׭^��ƭ�
		}

		sum_fitness = 0;
#pragma omp parallel for
		for (int individual = 0; individual < population_size; individual++)
		{
			sum_fitness += fitness[individual];
		}
		mean_fitness = sum_fitness / (float)population_size;
	}

		        
		if(goodgeneration == 50)
		{	
#pragma omp parallel
	{
			srand((unsigned int)time(NULL)^omp_get_thread_num());
#pragma omp for
			for(int i = 0;i < 3;i++)
			{
				individual = (int)floorf((float)(population_size - 1) * generate_rand_zero_to_one());

				fitness[individual] = keepfitness[i];
				for(gene = 0; gene < number_of_gene; gene++)
					unknowns[individual][gene] = keepunknowns[i][gene];
			}
	}
			goodgeneration = 0;
		}
}

void scaling(float constant)//�u���ഫ
{
	float a = 0, b = 0;

	if (constant != 0 && (fitness[fittest_individual] - mean_fitness) < 0)
	{
		a = (constant - 1)*mean_fitness / (fitness[fittest_individual] - mean_fitness);
		b = (1 - a)*mean_fitness;

		sum_fitness = 0;
#pragma omp parallel for
		for (individual = 0; individual < population_size; individual++)
		{
			fitness[individual] = a*fitness[individual] + b;
			sum_fitness += fitness[individual];
		}
		mean_fitness = sum_fitness / (float)population_size;
	}
}


void pre_process()
{
	reset();//�ܼƪ���l��
	initial_population();//��l�Ƹs��
	calculate_fitness();//�p��A����
	statistics();//�׭^���N
	scaling(scaling_constant);//�u���ഫ
}


int selection(int mate, float sum_fitness_temp)//���L�����
{
	float sum = 0, roulette = 0;
	individual = -1;
	roulette = generate_rand_zero_to_one()*sum_fitness_temp;
	do
	{
		individual++;
		sum += fitness[individual];
	} while (!((sum >= roulette) || (individual == population_size - 1)));
	mate = individual;
	return mate;
}

void FFD()
{
	float FFD_SD[number_of_gene];

		for(int f = 0;f < 4;f++)
		{
			FFD_fitness[f] = 0;
				if(FFDsample[gene][f] == 0)
				{
					BL_time = (int)(90 * unknowns[mate1][0] + 31);
					BL_pre = (int)(80 * unknowns[mate1][1] + 6);
					Decline = (int)(3 * unknowns[mate1][2] + 2);
					continued = (int)(5 * unknowns[mate1][3] + 2);
					
#pragma omp parallel
	{
		BLwm = (float*)malloc(sizeof(float) * BL_time);
		BLarr = (float*)malloc(sizeof(float) * (ROW-BL_time));
	#pragma omp for
					for(int p=0;p<datas;p++) {	
							for(int i=0; i!=BL_time; ++i) {
								BLwm[i] = 0;
							}
							for(int i=0; i!=(ROW-BL_time); ++i) {
								BLarr[i] = 0;
							}
								times=0,ODI=0,Judge=1;


						BL_top = (int)ceil(BL_time*((float)BL_pre/100));

						for(int e = BL_time ; e < ROW2[p] ; e++) {

							for(int m=e-BL_time ; m < e ; m++) {
								BLwm[e - m] = arri[p][m];
							}     
							sort(BLwm, BLwm+BL_time);
				
							top_mean = 0;
							for(int k=0 ; k < BL_top ; k++) { 
								top_mean += (BLwm[k + ((BL_time - 1) - BL_top)]);         
							}
							BLarr[e-BL_time] = top_mean / (float)BL_top;
						}


						for(int i=BL_time;i<ROW2[p];i++) {
							if((((arri[p][i]/BLarr[i-BL_time])-1.0000)*(-100))>=Decline)
							{
								times++;
							}
							else{ 
								times=0;
								Judge=1;
							}
							if((times>=continued))
							{
								if(Judge==1)
								{
									ODI++;
									Judge=0;
								}
							}
						}

						if(ROW2[p] != 0) {
							ODI_L[p]=(float)ODI/((float)ROW2[p]/3600);
						} else {
							ODI_L[p]=0;
						}
					}
	}
						xsum = 0, ysum = 0 , xysum = 0, x2sum = 0, y2sum = 0;
					
						for(int i=0;i< datas;i++) {
							xsum += ODI_L[i];
						}
						for(int i=0;i< datas;i++) {
							ysum += AHI[i];
						}

						xsum = xsum / ((float)datas);
						ysum = ysum / ((float)datas);

						for(int i=1;i< (int)datas;i++) {
								xysum += (ODI_L[i] - xsum) * (AHI[i] - ysum);
								x2sum += pow((ODI_L[i] - xsum), 2);
								y2sum += pow((AHI[i] - ysum), 2);
						}
							FFD_fitness[f] = (xysum / (sqrt(x2sum) * sqrt(y2sum)));

				}
				else
				{
					FFD_fitness[f] = 0;

					BL_time = (int)(90 * unknowns[mate2][0] + 31);
					BL_pre = (int)(80 * unknowns[mate2][1] + 6);
					Decline = (int)(3 * unknowns[mate2][2] + 2);
					continued = (int)(5 * unknowns[mate2][3] + 2);

					for(int p=0;p<datas;p++) {		
							for(int i=0; i!=BL_time; ++i) {
								BLwm[i] = 0;
							}
							for(int i=0; i!=(ROW-BL_time); ++i) {
								BLarr[i] = 0;
							}
								times=0,ODI=0,Judge=1;


						BL_top = (int)ceil(BL_time*((float)BL_pre/100));

						for(int e = BL_time ; e < ROW2[p] ; e++) {

							for(int m=e-BL_time ; m < e ; m++) {
								BLwm[e - m] = arri[p][m];
							}     
							sort(BLwm, BLwm+BL_time);
				
							top_mean = 0;
							for(int k=0 ; k < BL_top ; k++) { 
								top_mean += (BLwm[k + ((BL_time - 1) - BL_top)]);         
							}
							BLarr[e-BL_time] = top_mean / (float)BL_top;
						}


						for(int i=BL_time;i<ROW2[p];i++) {
							if((((arri[p][i]/BLarr[i-BL_time])-1.0000)*(-100))>=Decline)
							{
								times++;
							}
							else{ 
								times=0;
								Judge=1;
							}
							if((times>=continued))
							{
								if(Judge==1)
								{
									ODI++;
									Judge=0;
								}
							}
						}

						if(ROW2[p] != 0) {
							ODI_L[p]=(float)ODI/((float)ROW2[p]/3600);
						} else {
							ODI_L[p]=0;
						}
					}

						xsum = 0, ysum = 0 , xysum = 0, x2sum = 0, y2sum = 0;
					
						for(int i=0;i< datas;i++) {
							xsum += ODI_L[i];
						}
						for(int i=0;i< datas;i++) {
							ysum += AHI[i];
						}

						xsum = xsum / ((float)datas);
						ysum = ysum / ((float)datas);

						for(int i=1;i< (int)datas;i++) {
								xysum += (ODI_L[i] - xsum) * (AHI[i] - ysum);
								x2sum += pow((ODI_L[i] - xsum), 2);
								y2sum += pow((AHI[i] - ysum), 2);
						}
							FFD_fitness[f] = (xysum / (sqrt(x2sum) * sqrt(y2sum)));

				}
		}

#pragma omp parallel for
	for (gene = 0; gene < number_of_gene; gene++)
	{
		FFD_T = 0;
		FFD_F = 0;
		for (int f = 0; f < 32; f++)
		{
			FFD_T += FFD_fitness[f];
		}
		FFD_SD[gene] = abs(FFD_F - FFD_F);
		if (FFD_F < FFD_T)//��ܰ^�m�ȸ��j
		{
			new_unknowns[new_individual][gene] = unknowns[mate2][gene];
		}
		else
		{
			new_unknowns[new_individual][gene] = unknowns[mate1][gene];
		}

	}

	max_fitness = 10000000.0f;//SD�ȳ̤p

	for (gene = 0; gene < number_of_gene; gene++)
	{
		if (FFD_SD[gene] < max_fitness)
		{
			max_fitness = FFD_SD[gene];
			fittest_individual = gene;
		}
	}

#pragma omp parallel for
	for (gene = 0; gene < number_of_gene; gene++)//�̨�,���Φ]�l
	{
		new_unknowns[new_individual+1][gene] = new_unknowns[new_individual][gene];
	}

	if (new_unknowns[new_individual][fittest_individual] == unknowns[mate1][fittest_individual])
	{
		new_unknowns[new_individual + 1][fittest_individual] = unknowns[mate2][fittest_individual];
	}
	else
	{
		new_unknowns[new_individual + 1][fittest_individual] = unknowns[mate1][fittest_individual];
	}
}

void BLX()
{
	/*printf("BLX��...\n");*/
#pragma omp parallel
	{
		b_max = 0, b_min = 0, I = 0;
#pragma omp for
		for (gene = 0; gene < number_of_gene; gene++)
		{
			b_max = max(unknowns[mate1][gene], unknowns[mate2][gene]);
			b_min = min(unknowns[mate1][gene], unknowns[mate2][gene]);
			I = b_max - b_min;

			do
			{
				new_unknowns[new_individual][gene] = b_max - (float)generate_rand_zero_to_one()*(1 + 2 * alpha)*I;
			} while (!((new_unknowns[new_individual][gene] <= RANGE_UPPER) && (new_unknowns[new_individual][gene] >= RANGE_LOWER)));

			do
			{
				new_unknowns[new_individual + 1][gene] = b_min + (float)generate_rand_zero_to_one()*(1 + 2 * alpha)*I;
			} while (!((new_unknowns[new_individual + 1][gene] <= RANGE_UPPER) && (new_unknowns[new_individual + 1][gene] >= RANGE_LOWER)));
		}
	}
}

void cross_over()//�i���t
{
		for(new_individual = 0;new_individual < population_size;new_individual += 2)//!		/*Loop over the population choosing pairs of mates.*/
		{
			mate1 = selection(mate1, sum_fitness);//���L�����
			mate2 = selection(mate2, sum_fitness);//���L�����
			if (generate_rand_zero_to_one() <= cross_over_probability)//�����t
			{
				FFD();
			}
			else//�S����t
			{
#pragma omp parallel for
				for (gene = 0; gene < number_of_gene; gene++)
				{
					new_unknowns[new_individual][gene] = unknowns[mate1][gene];
					new_unknowns[new_individual + 1][gene] = unknowns[mate2][gene];
				}
			}
		}
}


void RM()
{
	float delta = 0.0f;
	
#pragma omp parallel
	{
		srand((unsigned int)time(NULL)^omp_get_thread_num());
#pragma omp for
		for (int individual = 0; individual < population_size; individual++)
		{
			for (gene = 0; gene < number_of_gene; gene++)
			{
				if (generate_rand_zero_to_one() < mutation_probability)
				{
					do
					{
						mutate_individual = 0.0f;
						delta = max(2 * (new_unknowns[individual][gene] - RANGE_LOWER), 2 * (RANGE_UPPER - new_unknowns[individual][gene]));
					} while (!((mutate_individual <= RANGE_UPPER) && (mutate_individual >= RANGE_LOWER)));

					new_unknowns[individual][gene] = mutate_individual;
				}
			}
		}
	}
}


void mutation()//����
{
	RM();
}



void replace()//���N
{
	#pragma omp parallel for
	for (individual = 0; individual < population_size; individual++)
	{
		for (gene = 0; gene < number_of_gene; gene++)
		{
			unknowns[individual][gene] = new_unknowns[individual][gene];
		}
	}
	calculate_fitness();//�p��A����
	statistics();//�׭^���N
	scaling(scaling_constant);//�u���ഫ
}

void ga_algorithm()//�i��GA�t��
{
	for (generation = 0; generation < max_generation; generation++)
	{
		
//�h�׭^�O�d
			if(goodgeneration == 49)
			{
				float lotMaxFitness = 100000000000.0f;
				for(int i = 0;i < 3;i++)
				{   
					max_fitness = -10000000000.0f;
					for(gene = 0; gene < number_of_gene; gene++)
					{
						if(fitness[individual] > max_fitness && fitness[individual] < lotMaxFitness)    /*�p�G��X���餤�̨θѪ̫h���N���̤j�̨θ�*/
						{
							max_fitness = fitness[individual];
							fittest_individual = individual;        /*�O���̨δ��խ���*/
						}
					}
					lotMaxFitness = max_fitness;
					keepfitness[i] = fitness[fittest_individual];
					for(gene = 0; gene < number_of_gene; gene++)
						keepunknowns[i][gene] = unknowns[fittest_individual][gene];
				}
			}

		cross_over();//�i���t
		mutation();//����
		replace();//���N

			if(generation % 50 == 0)
			{
				

				if(fopen_s(&fp,"ODI_L.txt","a") ==0){ /* open file pointer */

				BL_time = (int)((90 * elite_unknowns[0] + 31));
				BL_pre = (int)((80 * elite_unknowns[1] + 6));
				Decline = (int)((3 * elite_unknowns[2] + 2));
				continued = (int)((5 * elite_unknowns[3] + 2));

					for(gene = 0; gene < number_of_gene; gene++)
						fprintf(fp,"%5.8f ",elite_unknowns[gene]);
					fprintf(fp,"\n%5.8f\n",elite_fitness);
					fprintf(fp,"BL_time %i ",BL_time);
					fprintf(fp,"BL_pre %i ",BL_pre);
					fprintf(fp,"Decline %i ",Decline);
					fprintf(fp,"continued %i\n\n",continued);
		
			
#pragma omp parallel
	{
			BLwm = (float*)malloc(sizeof(float) * BL_time);
			BLarr = (float*)malloc(sizeof(float) * (ROW-BL_time));
		#pragma omp for
					for(int p=0;p<datas;p++) {		
							for(int i=0; i!=BL_time; ++i) {
								BLwm[i] = 0;
							}
							for(int i=0; i!=(ROW-BL_time); ++i) {
								BLarr[i] = 0;
							}
								times=0,ODI=0,Judge=1;


						BL_top = (int)ceil(BL_time*((float)BL_pre/100));

						for(int e = BL_time ; e < ROW2[p] ; e++) {

							for(int m=e-BL_time ; m < e ; m++) {
								BLwm[e - m] = arri[p][m];
							}     
							sort(BLwm, BLwm+BL_time);
				
							top_mean = 0;
							for(int k=0 ; k < BL_top ; k++) { 
								top_mean += (BLwm[k + ((BL_time - 1) - BL_top)]);         
							}
							BLarr[e-BL_time] = top_mean / (float)BL_top;
						}


						for(int i=BL_time;i<ROW2[p];i++) {
							if((((arri[p][i]/BLarr[i-BL_time])-1.0000)*(-100))>=Decline)
							{
								times++;
							}
							else{ 
								times=0;
								Judge=1;
							}
							if((times>=continued))
							{
								if(Judge==1)
								{
									ODI++;
									Judge=0;
								}
							}
						}

						if(ROW2[p]!= 0) {
							fprintf(fp,"%5.8f\n",(float)ODI/((float)ROW2[p]/3600));
						} else {
							fprintf(fp,"%i\n",0);
						}
					}
		}
					xsum = 0, ysum = 0 , xysum = 0, x2sum = 0, y2sum = 0;
				
					for(int i=0;i< datas;i++) {
						xsum += ODI_L[i];
					}
					for(int i=0;i< datas;i++) {
						ysum += AHI[i];
					}

					xsum = xsum / ((float)datas);
					ysum = ysum / ((float)datas);

					for(int i=1;i< (int)datas;i++) {
							xysum += (ODI_L[i] - xsum) * (AHI[i] - ysum);
							x2sum += pow((ODI_L[i] - xsum), 2);
							y2sum += pow((AHI[i] - ysum), 2);
					}
				printf(" %5.8f \n",(xysum / (sqrt(x2sum) * sqrt(y2sum))));
					fclose(fp);
				}
				
				for(gene = 0; gene < number_of_gene; gene++)
					printf("EliteUnknowns %5.8f \n",elite_unknowns[gene]);
				printf(" %5.8f \n",elite_fitness);
				printf("BL_time %i BL_pre %i Decline %i continued %i\n\n",BL_time,BL_pre,Decline,continued);
			}
			
//�����|�N
		if(restgenerate == terminate_generation)
		{
			generation = max_generation;
		}

	}
}


int main()
{
		printf("OSA���Ū����...\n");
		OSA_Load();
		printf("���Ū������\n");
		printf("OSA�̨Τƶ}�l����...\n");
		printf("�s��� %i   �|�N�� %i   ��t���v %5.3f  ���ܾ��v %5.3f \n\n",population_size,max_generation,cross_over_probability,mutation_probability);

	/*
#pragma omp parallel
	{
		srand((unsigned int)time(NULL)^omp_get_thread_num());*/
//#pragma omp for
			for (int i = 0; i < count_times; i++)
			{
				mean_fitness = 0;
				pre_process();//�e�B�z:�q�n��l��
				ga_algorithm();//�i��GA�t��
			}
	//}
	
	system("pause");
	return 0;
}