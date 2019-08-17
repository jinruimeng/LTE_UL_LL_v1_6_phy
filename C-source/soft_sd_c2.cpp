/* 
   This function implements a soft output sphere decoder.
   Based on the paper: C. Studer, M. Wenk, A. Burg, and H. Bölcskei: "Soft-Output Sphere
   Decoding: Performance and Implementation Aspects", Asilomar 2006
 
   R ... is an upper triangular matrix obtained from the QR decomposition of the MIMO channel H
   s ... received symbol (column) vector, s=Q^H*y 
   dist_ZF ... Distance for the zero forcing solution
   symbols_ZF ... ZF solution
   M ... number of bits encoded in every layer [1 x nTA]
   symbol_alphabet ... for the demapping [nTA x 2^max(M)], filled with zeros for smaller symbol alphabets
   bittable ... matrix containing the bits according to the symbol_alphabet [sum(M) x 2^max(M)]

   created: 10. April 2007, Christian Mehlführer, chmehl@nt.tuwien.ac.at
 */


static double inf = 99999999.9; // infinity

#define VISUAL_C 0

#if VISUAL_C==1
#include <stdlib.h>
#include <malloc.h>
#else
#include <matrix.h>
#endif

#if VISUAL_C==0
void soft_sd_c2( double *R_re, double *R_im,
				double *s_re, double *s_im, 
				double *dist_ZF,
				int *symbols_ZF_i,
				double *symbol_alphabet_re, double *symbol_alphabet_im,
				bool *bittable,
				double *LLR,
				int nSym, int *M, int nT, int nR, int total_bits) 
{
#else
void main(void)
{
	/*	double R_re[] = {-2.2076, 0.0,   0.4758, 2.1130};
	double R_im[] = {      0.0, 0.0, - 0.2040, 0.0};
	double s_re[] = {-1.7620, -1.4973};
	double s_im[] = {2.0327, 1.4896};
	double dist_ZF[] =  {1.8778e-004};
	int symbols_ZF_i[] = {2,  3};
	double symbol_alphabet_re[] = {0.7071,   0.7071,  -0.7071,  -0.7071,  0.7071,   0.7071,  -0.7071,  -0.7071};
	double symbol_alphabet_im[] = {0.7071,  -0.7071,   0.7071,  -0.7071,  0.7071,  -0.7071,   0.7071,  -0.7071};
	bool bittable[] = {0,0,0,0,1,0,1,0,0,1,0,1,1,1,1,1}; 

	int M[] = {2,2};
	double LLR[4];
	int nSym=1;
	int nT=2;
	int nR=2;
	int total_bits = 4; */

	double R_re[] = {-3.5085, 0, 0, -0.3170, 0.6803, 0, 0.7190, 0.7202, 0.7763};
	double R_im[] = {0, 0, 0, 0.1131, 0, 0, -0.7122, 0.5091, 0};
	double s_re[] = {-4.7642, -3.8677, 7.2334};
	double s_im[] = {1.4122, -16.3026, -3.4942};
	double dist_ZF[] =  {362.9078};
	int symbols_ZF_i[] = {2, 2, 1};
	double symbol_alphabet_re[] = {0.7071, 0.7071, -0.7071, -0.7071, 1.0000, -1.0000, 0, 0, 1.0000, -1.0000, 0, 0};
	double symbol_alphabet_im[] = {0.7071, -0.7071, 0.7071, -0.7071, 0, 0, 0, 0, 0, 0, 0, 0};
	bool bittable[] = {0,0,0,0,1,0,1,1,0,1,0,0,1,1,0,0}; 

	int M[] = {2,1,1};
	double LLR[4];
	int nSym=1;
	int nT=3;
	int nR=3;
	int total_bits = 4; 
#endif

	int max_layer, max_Nr_of_symbols, current_layer, kk, mm, ss, *current_symbol_i, *ML_symbol_i, index, index1, sym, M_sum, *Nr_of_symbols; 
	double lambda_ML, *lambda_CH, *dist, *s_re_cur, *s_im_cur, dist_temp_re, dist_temp_im, lambda_CH_max;
	bool searching_criterion, end_while;

	max_layer = nT;              /* number of layers (transmit, receive antennas) */

	/* allocate memory */
#if VISUAL_C==1
	lambda_CH        = (double *)malloc( total_bits*sizeof(double) ); /* distances of the counter hypotheses */
	current_symbol_i = (int *)malloc( max_layer*sizeof(int) );        /* current symbol values: 1...2^M */
	ML_symbol_i      = (int *)malloc( max_layer*sizeof(int) );        /* ML solution 1...2^M, initialize to ZF solution */
	dist             = (double *)malloc( max_layer*sizeof(double) );  /* current distances for all layers */
	Nr_of_symbols    = (int *)malloc( max_layer*sizeof(int) );        /* number of symbols in every layer */
#else
	lambda_CH     = (double *)mxCalloc( total_bits, sizeof(double) ); /* distances of the counter hypotheses */
	current_symbol_i = (int *)mxCalloc( max_layer, sizeof(int) );     /* current symbol values: 1...2^M */
	ML_symbol_i = (int *)mxCalloc( max_layer, sizeof(int) );          /* ML solution 1...2^M, initialize to ZF solution */
	dist = (double *)mxCalloc( max_layer, sizeof(double) );           /* current distances for all layers */
	Nr_of_symbols = (int *)mxCalloc( max_layer, sizeof(int) );        /* number of symbols in every layer */
#endif

	max_Nr_of_symbols = 0;
	M_sum = 0;
	for(kk=0; kk<max_layer; kk++)
	{
		Nr_of_symbols[kk] = 2;
		for(mm=1; mm<M[kk]; mm++)
			Nr_of_symbols[kk] *= 2;
		M_sum += M[kk];
		if (Nr_of_symbols[kk]>max_Nr_of_symbols)
			max_Nr_of_symbols = Nr_of_symbols[kk];
	}

	for(ss=0; ss<nSym; ss++)
	{
		/* initialize variables */
		lambda_ML = dist_ZF[ss];            /* distance of the ML solution, initialize to ZF solution */
		current_layer = max_layer;          /* current layer */
		for(kk=0; kk<total_bits; kk++) 
			lambda_CH[kk] = inf;
		for(kk=0; kk<max_layer; kk++)
		{
			current_symbol_i[kk] = 1;
			ML_symbol_i[kk] = symbols_ZF_i[kk+ss*max_layer];
			dist[kk] = 0;
		}
		s_re_cur = &s_re[ss*max_layer];		/* set pointer to received symbol vector */
		s_im_cur = &s_im[ss*max_layer];		/* set pointer to received symbol vector */
		searching_criterion = true;

		/* start the number crunching */
		while (searching_criterion)
		{
			index = current_layer-1+(max_layer-1)*max_layer;
			sym = max_Nr_of_symbols*(max_layer-1)+current_symbol_i[max_layer-1]-1;
			dist_temp_re = s_re_cur[current_layer-1] - (R_re[index]*symbol_alphabet_re[sym] - R_im[index]*symbol_alphabet_im[sym]);
			dist_temp_im = s_im_cur[current_layer-1] - (R_im[index]*symbol_alphabet_re[sym] + R_re[index]*symbol_alphabet_im[sym]);

			for (kk=max_layer-1; kk>=current_layer; kk--)
			{
				index = current_layer-1+(kk-1)*max_layer;
				sym = max_Nr_of_symbols*(kk-1)+current_symbol_i[kk-1]-1;
				dist_temp_re -= R_re[index]*symbol_alphabet_re[sym] - R_im[index]*symbol_alphabet_im[sym];
				dist_temp_im -= R_im[index]*symbol_alphabet_re[sym] + R_re[index]*symbol_alphabet_im[sym];
			}

			if (current_layer==max_layer)
				dist[current_layer-1] = dist_temp_re*dist_temp_re + dist_temp_im*dist_temp_im;
			else
				dist[current_layer-1] = dist[current_layer] + dist_temp_re*dist_temp_re + dist_temp_im*dist_temp_im;

			if (current_layer > 1)
			{
				if (dist[current_layer-1] < lambda_ML)
					current_layer--;	/* proceed to next layer */
				else
				{
					lambda_CH_max = 0;
					index =0;
					for (kk=0; kk<max_layer; kk++)      /* pruning condition */
					{
						for (mm=0; mm<M[kk]; mm++)
						{
							if (kk<current_layer-1)
							{
								if (lambda_CH[index]>lambda_CH_max)
									lambda_CH_max = lambda_CH[index];
							}
							else
								if ((lambda_CH[index]>lambda_CH_max) && (bittable[index+M_sum*(ML_symbol_i[kk]-1)]!=bittable[index+M_sum*(current_symbol_i[kk]-1)]))
									lambda_CH_max = lambda_CH[index];
							index++;
						}
					}
					if (dist[current_layer-1] < lambda_CH_max)
						current_layer--;	/* proceed to next layer */
					else
					{
						end_while = false;
						while (!end_while)
							if (current_symbol_i[current_layer-1] == Nr_of_symbols[current_layer-1])
								if (current_layer < max_layer)
								{
									current_symbol_i[current_layer-1] = 1;      /* start with first symbol */
									current_layer++;							/* in the next layer */
								}
								else
								{
									searching_criterion = false;				/* terminate search */
									end_while = true;
								}
							else
							{
								current_symbol_i[current_layer-1]++;				/* go to next symbol */
								end_while = true;
							}	
					}
				}
			}        

			else        /* bottom layer reached (current_layer==1) */
				if (dist[current_layer-1] < lambda_ML)		/* new hypothesis found */
				{
					index = 0;
					for (kk=0; kk<max_layer; kk++)			/* update counter hypotheses distances */
						for (mm=0; mm<M[kk]; mm++)
						{
							if (bittable[index+M_sum*(ML_symbol_i[kk]-1)]!=bittable[index+M_sum*(current_symbol_i[kk]-1)])
								lambda_CH[index] = lambda_ML;
							index++;
						}
						for (kk=0; kk<max_layer; kk++)			/* update ML symbol */
							ML_symbol_i[kk] = current_symbol_i[kk];
						lambda_ML = dist[current_layer-1];		/* overwrite ML distance */
				}
				else										/* check counter hypotheses distances */
				{
					index = 0;
					for (kk=0; kk<max_layer; kk++)			/* update counter hypotheses distances */
						for (mm=0; mm<M[kk]; mm++)
						{
							if ((lambda_CH[index]>dist[current_layer-1])&&(bittable[index+M_sum*(ML_symbol_i[kk]-1)]!=bittable[index+M_sum*(current_symbol_i[kk]-1)]))
								lambda_CH[index] = dist[current_layer-1];
							index++;
						}
						end_while = false;
						while (!end_while)
							if (current_symbol_i[current_layer-1] == Nr_of_symbols[current_layer-1])
								if (current_layer < max_layer)
								{
									current_symbol_i[current_layer-1] = 1;      /* start with first symbol */
									current_layer++;							/* in the next layer */
								}
								else
								{
									searching_criterion = false;				/* terminate search */
									end_while = true;
								}
							else
							{
								current_symbol_i[current_layer-1]++;				/* go to next symbol */
								end_while = true;
							}
				}
		}

		/* calculate LLR values from distances  */
		index1 = ss*total_bits;
		index = 0;
		for (kk=0; kk<max_layer; kk++)
		{
			for (mm=0; mm<M[kk]; mm++)
			{
				if (!bittable[index+M_sum*(ML_symbol_i[kk]-1)])
					LLR[index1] = lambda_ML - lambda_CH[index];
				else
					LLR[index1] = lambda_CH[index] - lambda_ML;
				index++;
				index1++;
			}
		}
	}
	/* free memory    */
#if VISUAL_C==1
	free(lambda_CH);
	free(current_symbol_i);
	free(ML_symbol_i);
	free(dist);
#else
	mxFree(lambda_CH);
	mxFree(current_symbol_i);
	mxFree(ML_symbol_i);
	mxFree(dist);
#endif

	/*	mexPrintf("LLR[0]=%f, LLR[1]=%f, LLR[2]=%f, LLR[3]=%f\n",  */
	/*		LLR[0], LLR[1], LLR[2], LLR[3]); */
} 
