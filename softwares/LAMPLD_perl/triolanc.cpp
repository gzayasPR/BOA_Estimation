#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include "hmm_phase.h"
#include "random.h"
#include <vector>
#include <iostream>
using namespace std ;


//#include <gsl/gsl_statistics_double.h>

double totalProbUno_forward_fast(char* ch, HAP_HMM hmmOne, HAP_HMM hmmTwo,int start_pos);
double totalProbUno_backward_fast(char* ch, HAP_HMM hmmOne, HAP_HMM hmmTwo,int start_pos);
/*variables*/
double *pi_temp;
double*** trans_temp;
STATE**	states_temp;
double **forward_hap_array;
double **backward_hap_array;
double **array_uno_precomp1;
double *** totalProbUnoMatrix_forward;
double *** totalProbUnoMatrix_backward;
double ** totalProbUnoMatrix_unk;
double *** totalProbMatrix; // Genotype Total matrix for single genotype

ANC_HMM *anc_map;
char **admx_samples;
int nr_admx_samples;
char **admx_samples_anc;
char **admx_samples_inf_anc;
char **admx_samples_inf_phase;

double ****admx_inf_array;


double eps = 0.1;
double founder_bias = .9;
int max_steps = 100;
int hmmCount = 3;
double noise = .3;
int *true_sol;
int *inf_sol;
int window_padding = 0;
int nr_anc = 2;
int nr_total_snps;
int *snps_pos;
double *snps_map;
int num_gen;
double recomb_rate;

/* defs*/
#define MAX_SNPS 49600
#define MAJOR_EM_PROB 0.9
#define ADD_NOISE 1
#define LOGZERO -100000 
#define MAXBLOCK 600
//10 maximum switches
#define MAXNRPHASE 1024
#define MAXWINS 20000
#define NRANC 3
#define MAXAMBGPOS 2


typedef struct{
	char hap[4][MAXBLOCK];
	int ancestries[4];
	double prob;	
}TRIOPHASE;

int starting_snp = 0;

/*functions*/
void fill_forward_backward_hap_arrays(char* hap ,HAP_HMM hmm,int start_pos);
void read_training_haps(int ancestry, const char* filename);
void read_samples(const char* filename);
void read_true_sol(const char* filename);
void read_pos_file(const char* filename);
void read_map_file(const char* filename);
void trainHMM(int h,int start_pos);
void alloc_mem_other();
void get_tot_acc(int);
int convAlleleCharToInteger(char val);
double mylog(double x){if(x<=0) return 0; return  log(x);};
void fix_ancestry(int right_snp);
void fix_ancestry_hap(int right_snp);
char* getgeno(char* hap1, char* hap2, int nr_snp);
double gethapprob(char* hap ,HAP_HMM hmm, int start_pos);

inline double eexp(double x){if(x==LOGZERO) return 0; return exp(x);}
inline double eln(double x){ if(x==0) return LOGZERO; return log(x);}
inline double elnsum(double elnx, double elny){
	if(elnx==LOGZERO || elny ==LOGZERO){
		if(elnx==LOGZERO){
			return elny;
		}else{
			return elnx;
		}
	}else{
		if(elnx > elny){
			return elnx + eln(1+exp(elny-elnx));
		}else{
			return elny + eln(1+exp(elnx-elny));
		}
	}
} 

inline double elnproduct(double elnx, double elny){
	if(elnx==LOGZERO || elny ==LOGZERO){
		return LOGZERO;
	}else{
		return elnx+elny;  
	}
}
double getTransitionProb(int prev_a1, int a1, 
						 int prev_a2, int a2,
						 int prev_a3, int a3,
						 int prev_a4, int a4,
						 double genDist, int numGen, double *alpha){
	
	double prob = 0;
	if(a1!= prev_a1){
		prob += log( (numGen - 1)*genDist*alpha[a1] );
	}
	if(a2!= prev_a2){
		prob += log( (numGen - 1)*genDist*alpha[a2] );
	}
	if(a3!= prev_a3){
		prob += log( (numGen - 1)*genDist*alpha[a3] );
	}
	if(a4!= prev_a4){
		prob += log( (numGen - 1)*genDist*alpha[a4] );
	}
	
	return (prob);	
}

int main(int argc, char** argv){
	
  double alpha[3];
  int tot_c,tot_e;
  double anc_gen_forward[3][3][MAX_SNPS];
  double anc_gen_backward[3][3][MAX_SNPS];
  double bkpt_prob[MAX_SNPS];
  bool isAmbg[MAX_SNPS];
  int best_s = -1;
  double best_l;
  int left[2],right[2],prev_win_r[2];
  int best_no_bkpt[2];
  int start_pos = 0;
  char haps[4][MAXBLOCK];
  int numPhases[MAXWINS];
  int numAmbg[MAXWINS];
  int startingSnp[MAXWINS];
  int endingSnp[MAXWINS];
  float ******lik_array;//[MAXWINS][MAXNRPHASE][NRANC][NRANC][NRANC][NRANC];
  float ******forward_array;//[MAXWINS][MAXNRPHASE][NRANC][NRANC][NRANC][NRANC];
  int viterbi_phase[MAXWINS];
  int viterbi_anc[MAXWINS][4];
  int x[20];
  int lastwin =0;
  int winsize;
  int increment;
  int maxAmbgWindowPos = MAXAMBGPOS;
  nr_founders = 8;
  hmmCount = 3;
  alpha[0]=alpha[1]=alpha[2] = 0.33;
  nr_anc = 3;
  recomb_rate = 0.00000001;
  num_gen = 7;
  
  init_rand(1);
  anc_map  = (ANC_HMM*)calloc(nr_anc, sizeof(ANC_HMM));
  
  
  //read_true_sol("haps.anc");
  //read_map_file("/home/bp45/projects/localanc/latino-trio/latino/chr1.genmap");
  
  nr_snp = 50;
  window_padding = 0;
  increment = nr_snp; 
  
  
  if(argc != 10){
    fprintf(stderr,"USAGE:: %s <ambg-pos-per-win(<=10)> <nr-founders> <pos-file> <pop0-haps> <pop1-haps> <pop2-haps> <trio-file(m/f/c)> <file-out> <smoothing/0-no,1-yes>\n", argv[0]);
    exit(0);
  }
  bool smoothing = false;
  if(atol(argv[9]) > 0){
    smoothing = true;
  }
  maxAmbgWindowPos  = atol(argv[1]);
  nr_founders = atol(argv[2]);
  read_pos_file(argv[3]);
  read_training_haps(0,argv[4]);
  read_training_haps(1,argv[5]);
  read_training_haps(2,argv[6]);
  read_samples(argv[7]);


  
  fprintf(stderr,"Read %d anc0 %d anc1 %d anc2 training haplotypes\n",anc_map[0].nr_hap,anc_map[1].nr_hap,anc_map[2].nr_hap);
  fprintf(stderr,"Read %d admixed samples\n",nr_admx_samples);
  fprintf(stderr,"Using window size %d with padding %d\n",nr_snp, window_padding);
	
	
  lik_array = (float******) calloc(MAXWINS,sizeof(float*****));
  forward_array = (float******) calloc(MAXWINS,sizeof(float*****));
  for (int i=0; i < MAXWINS; i++){
    lik_array[i] = (float*****) calloc(MAXNRPHASE,sizeof(float****));
    forward_array[i] = (float*****) calloc(MAXNRPHASE,sizeof(float****));
    for(int a1=0;a1<NRANC;a1++){
      lik_array[i][a1] = (float****) calloc(NRANC,sizeof(float***));
      forward_array[i][a1] = (float****) calloc(NRANC,sizeof(float***));
      for(int a2=0;a2<NRANC;a2++){
	lik_array[i][a1][a2] = (float***) calloc(NRANC,sizeof(float**));
	forward_array[i][a1][a2] = (float***) calloc(NRANC,sizeof(float**));
	for(int a3=0;a3<NRANC;a3++){
	  lik_array[i][a1][a2][a3] = (float**) calloc(NRANC,sizeof(float*));
	  forward_array[i][a1][a2][a3] = (float**) calloc(NRANC,sizeof(float*));
	  for(int a4=0;a4<NRANC;a4++){
	    //lik_array[i][a1][a2][a3][a4] = (float*) calloc(NRANC,sizeof(float));
	    //lik_array[i][j][a1][a2][a3][a4] = 0;
	    //lik_array[i][a1][a2][a3][a4] = (float*) calloc (MAXNRPHASE,sizeof(float));
	    //forward_array[i][a1][a2][a3][a4] = (float*) calloc (MAXNRPHASE,sizeof(float));
	    //fprintf(stderr,"%d %d %d %d %d %d %f\n",i,j,a1,a2,a3,a4,lik_array[i][j][a1][a2][a3][a4]);
	  }
	}
      }
    }
  }
  
  //winsize = (int)(sqrt(eps)/(2*(num_gen-1)*recomb_rate*(1-alpha[0]*alpha[0])*(1-alpha[1]*alpha[1])));
  
  //fprintf(stderr,"Using window size %d bp\n",winsize);
  
  //FILE *flog = fopen("log.txt","w");
  
  for(int h = 0; h < hmmCount; h++){
    alloc_mem_hmm(&anc_map[h].hmm);
  }
	
  alloc_mem_other();
  
  for(int h = 0; h < hmmCount; h++){
    anc_map[h].forward_hap_array = (double**)malloc(MAXWINS*sizeof(double*));
    for(int i =0; i < MAXWINS; i++){
      anc_map[h].forward_hap_array[i] = (double*)malloc(nr_founders*sizeof(double));
    }	
    anc_map[h].backward_hap_array = (double**)malloc(MAXWINS*sizeof(double*));
    for(int i =0; i < MAXWINS; i++)    {
	    anc_map[h].backward_hap_array[i] = (double*)malloc(nr_founders*sizeof(double));
    }
  }
  
  for (int ii = starting_snp; ii < nr_total_snps; ii++) {
    //assume no missing data for now!!
    if( admx_samples[0][ii] == '?' ||
	admx_samples[1][ii] == '?' ||
	admx_samples[2][ii] == '?' ||
	admx_samples[3][ii] == '?' 
	){
      isAmbg[ii] =  true;
    }else{
      isAmbg[ii] =  false;
    }
  }
  tot_c=tot_e=0;
  start_pos=starting_snp;
  lastwin=0;
  int window_idx = 0;
  //nr_total_snps = 5000;
  //fprintf(stderr,"%s\n",anc_map[0].hap_array[1]
  while(start_pos  < nr_total_snps){
    
    int end_pos = start_pos;
    
    end_pos = start_pos;
    int nr_ambg = 0;

    while(nr_ambg < maxAmbgWindowPos && end_pos < nr_total_snps && end_pos-start_pos < MAXBLOCK-20){
      if(isAmbg[end_pos]){
	nr_ambg++;
      }
      end_pos++;
    }

    //end_pos+=3;
    
    //    fprintf(stderr,"%d %d %d\n",start_pos,end_pos,nr_ambg); 
    nr_snp = end_pos - start_pos;
    
	  
    if((end_pos) > nr_total_snps){
      nr_snp = nr_total_snps-start_pos;
      end_pos = nr_total_snps;
      lastwin=1;
    }
    fprintf(stdout,"Admixed calculations at SNP Win [%d %d] %d real[%d %d (%d)]\n",start_pos,start_pos + nr_snp,nr_snp,
	    snps_pos[start_pos],snps_pos[start_pos+nr_snp],snps_pos[start_pos+nr_snp]-snps_pos[start_pos]);
    
    for(int h = 0; h < hmmCount; h++){
      trainHMM(h,start_pos);
    }
    //fprintf(stderr,"training complete\n");
    int ambg_pos = 0;
    for (int ii = start_pos; ii < end_pos; ii++) {
      if( isAmbg[ii] ){
	ambg_pos++;
	for(int iii = 0; iii<4 ;iii++){
	  haps[iii][ii-start_pos] = '?';
	}
      }else{
	for(int iii = 0; iii<4 ;iii++){
	  haps[iii][ii-start_pos] = admx_samples[iii][ii];
	}
	//isAmbg[ii] =  false;
      }
    }
    //fprintf(stderr,"ambg_pos %d \n",ambg_pos);

    //ambg_pos = 0;
    numAmbg[window_idx] = ambg_pos;
    numPhases[window_idx] = (int)pow(2,ambg_pos);
    startingSnp[window_idx] = start_pos;
    endingSnp[window_idx] = end_pos-1;
    

    
    if(ambg_pos > 10){
      fprintf(stderr,"increase array size %d!!\n",ambg_pos);
      exit(0);
    }
    for(int a1 =0; a1 < nr_anc; a1++){
      for(int a2 =0; a2 < nr_anc; a2++){
	for(int a3 =0; a3 < nr_anc; a3++){
	  for(int a4 =0; a4 < nr_anc; a4++){
	    lik_array[window_idx][a1][a2][a3][a4] = (float*) calloc (numPhases[window_idx],sizeof(float));
	    forward_array[window_idx][a1][a2][a3][a4] = (float*) calloc (numPhases[window_idx],sizeof(float));
	  }
	}
      }
    }
    for(int tt = 0 ; tt < pow(2,ambg_pos); tt++){
      //fprintf(stderr,"%d\n",tt);
      int ttt = tt;
      for(int idx = ambg_pos-1; idx>=0;idx--){
	x[idx] = ttt%2;
	ttt = ttt/2;
      }
      //fprintf(stderr,"%d\n",tt);
      int ambg_idx = 0;
      for (int ii = start_pos; ii < end_pos; ii++) {
	if( isAmbg[ii]){
	  if(x[ambg_idx] == 0){
	    haps[0][ii-start_pos] = '0' ;
	    haps[1][ii-start_pos] = '1' ;
	    haps[2][ii-start_pos] = '1' ;
	    haps[3][ii-start_pos] = '0' ;
	  }else{
	    haps[0][ii-start_pos] = '1' ;
	    haps[1][ii-start_pos] = '0' ;
	    haps[2][ii-start_pos] = '0' ;
	    haps[3][ii-start_pos] = '1' ;
	  }
	  ambg_idx++;
	}
      }

      double xx[4][NRANC];
      for(int a1 =0; a1 < nr_anc; a1++){
	for(int i =0; i < 4; i++){
	  xx[i][a1] = gethapprob(haps[i],anc_map[a1].hmm,  0);
	}
      }
      
      
      for(int a1 =0; a1 < nr_anc; a1++){
	for(int a2 =0; a2 < nr_anc; a2++){
	  for(int a3 =0; a3 < nr_anc; a3++){
	    for(int a4 =0; a4 < nr_anc; a4++){
	      double x1 = xx[0][a1];
	      double x2 = xx[1][a2];
	      double x3 = xx[2][a3];
	      double x4 = xx[3][a4];
	      double p = log(x1)+log(x2)+log(x3)+log(x4);
	      lik_array[window_idx][a1][a2][a3][a4][tt] = p;
	      //fprintf (stderr, "%d %d (%d %d %d %d) %e %e %e %e %f\n",window_idx,tt,a1,a2,a3,a4,x1,x2,x3,x4,lik_array[window_idx][a1][a2][a3][a4][tt]);
	    }
	  }
	}
      }
    }
    if(window_idx > 10000){
      fprintf(stderr,"increase array size!!\n");
      exit(0);
    }
    window_idx++;
    start_pos = end_pos;		
    if(lastwin){
      start_pos = nr_total_snps+1;
    }		
  }

  int nr_windows = window_idx;
  fprintf(stderr,"%d\n",nr_windows);
  for(int i =0; i < nr_windows;i++){
    fprintf(stderr,"%d\n",i);
    if(i==0){
      for(int j = 0; j< numPhases[i];j++)
	for(int a1 =0; a1 < nr_anc; a1++)
	  for(int a2 =0; a2 < nr_anc; a2++)
	    for(int a3 =0; a3 < nr_anc; a3++)
	      for(int a4 =0; a4 < nr_anc; a4++)
		forward_array[i][a1][a2][a3][a4][j] = lik_array[i][a1][a2][a3][a4][j];
    }else{						
      for(int j = 0; j< numPhases[i];j++)
	for(int a1 =0; a1 < nr_anc; a1++)
	  for(int a2 =0; a2 < nr_anc; a2++)
	    for(int a3 =0; a3 < nr_anc; a3++)
	      for(int a4 =0; a4 < nr_anc; a4++){
		
		double best_p = -9999999;
		int best_j = 0;
		int best_a1 = 0;
		int best_a2 = 0;
		int best_a3 = 0;
		int best_a4 = 0;
		double genDist = ( snps_pos[ startingSnp[i]] - snps_pos[endingSnp[i-1]])*recomb_rate;
		double trans,prev_lik,cur_lik;
		cur_lik = lik_array[i][a1][a2][a3][a4][j];
		
		
		for(int prev_j = 0; prev_j< numPhases[i-1];prev_j++)
		  {
		    
		    for(int prev_a1 =0; prev_a1 < nr_anc; prev_a1++)
		      for(int prev_a2 =0; prev_a2 < nr_anc; prev_a2++)
			for(int prev_a3 =0; prev_a3 < nr_anc; prev_a3++)
			  for(int prev_a4 =0; prev_a4 < nr_anc; prev_a4++)
			    {
			      trans = getTransitionProb(prev_a1,a1,prev_a2,a2,prev_a3,a3,prev_a4,a4,
							genDist,num_gen,alpha);
			      prev_lik = forward_array[i-1][prev_a1][prev_a2][prev_a3][prev_a4][prev_j];
			      
			      /*fprintf(stderr,"%d (%d,%d,%d,%d,%d) -> (%d,%d,%d,%d,%d) %f %f %f\n",
				i,prev_j,prev_a1,prev_a2,prev_a3,prev_a4,
				j,a1,a2,a3,a4,trans,prev_lik,cur_lik);
			      */
			      if(trans+prev_lik+cur_lik > best_p){
				best_p = trans+prev_lik+cur_lik;
				best_a1 = prev_a1;
				best_a2 = prev_a2;
				best_a3 = prev_a3;
				best_a4 = prev_a4;
				best_j = prev_j;
			      }
			    }
		  }
		forward_array[i][a1][a2][a3][a4][j] = best_p;
		/*fprintf(stderr,"%d (%d,%d,%d,%d,%d) -> (%d,%d,%d,%d,%d) %f\n",
		  i,best_j,best_a1,best_a2,best_a3,best_a4,
		  j,a1,a2,a3,a4,best_p);
		*/
	      }
		}
  }
  /*get max and trace back*/
  
  double lik = -9999999;
  int cur_a1 = 0;
  int cur_a2 = 0;
  int cur_a3 = 0;
  int cur_a4 = 0;
  int cur_j = 0;
  
  for(int i = nr_windows-1;  i >= 0 ;i--){
    int best_j = 0;
    int best_a1 = 0;
    int best_a2 = 0;
    int best_a3 = 0;
    int best_a4 = 0;
    lik = -99999999;
    if(i==nr_windows-1){
      for(int j = 0; j< numPhases[i];j++)
	for(int a1 =0; a1 < nr_anc; a1++)
	  for(int a2 =0; a2 < nr_anc; a2++)
	    for(int a3 =0; a3 < nr_anc; a3++)
	      for(int a4 =0; a4 < nr_anc; a4++){
		if(forward_array[i][a1][a2][a3][a4][j] > lik){
		  lik = forward_array[i][a1][a2][a3][a4][j];
		  best_j = j;
		  best_a1 = a1;
		  best_a2 = a2;
		  best_a3 = a3;
		  best_a4 = a4;
		  //fprintf(stderr,"%d (%d,%d,%d,%d,%d) %f\n",i,best_j,best_a1,best_a2,best_a3,best_a4, lik);getchar();
		}
	      }
      cur_a1 = best_a1;
      cur_a2 = best_a2;
      cur_a3 = best_a3;
      cur_a4 = best_a4;
      cur_j = best_j;
    }else {
      
      lik = -9999999;
      double genDist = ( snps_pos[ startingSnp[i+1]] - snps_pos[endingSnp[i]])*recomb_rate;
      double trans,prev_lik,cur_lik;
      cur_lik = lik_array[i+1][cur_a1][cur_a2][cur_a3][cur_a4][cur_j];
			
      for(int prev_j = 0; prev_j< numPhases[i];prev_j++)
	for(int prev_a1 =0; prev_a1 < nr_anc; prev_a1++)
	  for(int prev_a2 =0; prev_a2 < nr_anc; prev_a2++)
	    for(int prev_a3 =0; prev_a3 < nr_anc; prev_a3++)
	      for(int prev_a4 =0; prev_a4 < nr_anc; prev_a4++){
		
		trans = getTransitionProb(prev_a1,cur_a1,prev_a2,cur_a2,prev_a3,cur_a3,prev_a4,cur_a4,
					  genDist,num_gen,alpha);
		prev_lik = forward_array[i][prev_a1][prev_a2][prev_a3][prev_a4][prev_j];
		
		/*fprintf(stderr,"%d (%d,%d,%d,%d,%d) -> (%d,%d,%d,%d,%d) %f %f %f\n",
		  i,prev_j,prev_a1,prev_a2,prev_a3,prev_a4,
		  cur_j,cur_a1,cur_a2,cur_a3,cur_a4,trans,prev_lik,cur_lik);
								 */
		if(trans+prev_lik+cur_lik > lik){
		  lik = trans+prev_lik+cur_lik;
		  best_a1 = prev_a1;
		  best_a2 = prev_a2;
		  best_a3 = prev_a3;
		  best_a4 = prev_a4;
		  best_j = prev_j;
		}
	      }
      cur_a1 = best_a1;
      cur_a2 = best_a2;
      cur_a3 = best_a3;
      cur_a4 = best_a4;
      cur_j = best_j;
      
    }
    viterbi_phase[i] = cur_j;
    viterbi_anc[i][0] = cur_a1;
    viterbi_anc[i][1] = cur_a2;
    viterbi_anc[i][2] = cur_a3;
    viterbi_anc[i][3] = cur_a4;		
  }
  
  for( int i = 0 ; i < nr_windows; i++){
		
    //fprintf(stderr,"%d (%d\t%d,%d,%d,%d)\n",i,viterbi_phase[i],viterbi_anc[i][0],viterbi_anc[i][1],viterbi_anc[i][2],viterbi_anc[i][3]);
    
    for (int ii = startingSnp[i];  ii <= endingSnp[i]; ii++) {
      for(int i_hap = 0; i_hap < 4; i_hap++){
	admx_samples_inf_anc[i_hap][ii] =(char) (viterbi_anc[i][i_hap]+'0');
      }								
    }		
    int ttt = viterbi_phase[i];
    for(int idx = numAmbg[i]-1; idx>=0;idx--){
      x[idx] = ttt%2;
      ttt = ttt/2;
    }
    int ambg_idx = 0;
    for (int ii = startingSnp[i];  ii <= endingSnp[i]; ii++) {
      if( isAmbg[ii]){
	admx_samples_inf_phase[0][ii] = '?' ;
	admx_samples_inf_phase[1][ii] = '?' ;
	admx_samples_inf_phase[2][ii] = '?' ;
	admx_samples_inf_phase[3][ii] = '?' ;
	if(x[ambg_idx] == 0){
	  admx_samples_inf_phase[0][ii] = '0' ;
	  admx_samples_inf_phase[1][ii] = '1' ;
	  admx_samples_inf_phase[2][ii] = '1' ;
	  admx_samples_inf_phase[3][ii] = '0' ;
	}else{
	  admx_samples_inf_phase[0][ii] = '1' ;
	  admx_samples_inf_phase[1][ii] = '0' ;
	  admx_samples_inf_phase[2][ii] = '0' ;
	  admx_samples_inf_phase[3][ii] = '1' ;
	}
	ambg_idx++;
      }else {
	for(int i_hap = 0; i_hap < 4; i_hap++){
	  admx_samples_inf_phase[i_hap][ii] = admx_samples[i_hap][ii];				
	}											
      }
    }
  }
  //test this before using it
  if(smoothing){
    int ambgPos=0;
    double forwardArray[2][200];
    double backwardArray[2][200];
    int anc[2];
    for(int i_hap = 0; i_hap < 4; i_hap++){
      int i_snp = 0;
      while( i_snp < nr_total_snps){
	if(i_snp > 50 && i_snp < nr_total_snps-50 && admx_samples_inf_anc[i_hap][i_snp] != admx_samples_inf_anc[i_hap][i_snp-1]){
	  nr_snp = 100;
	  fprintf(stderr,"%d %d (%c -> %c)\n",i_hap,i_snp,admx_samples_inf_anc[i_hap][i_snp-1],admx_samples_inf_anc[i_hap][i_snp] );
	  /*found a breakpoint -> readjusting it*/
	  anc[0] = (int)(admx_samples_inf_anc[i_hap][i_snp-1]-'0');
	  anc[1] = (int)(admx_samples_inf_anc[i_hap][i_snp]-'0');
	  //fprintf(stderr,"%d %d\n",anc[0],anc[1]);
	  //getchar();
	  for(int h = 0; h < hmmCount; h++){
	    trainHMM(h, i_snp -50);
	  }	  
	  for(int i_anc = 0 ; i_anc < 2; i_anc++){			 
	    fill_forward_backward_hap_arrays(admx_samples_inf_phase[i_hap], anc_map[ anc[i_anc] ].hmm,  i_snp-50);
	    for(int ii = 0; ii < 100; ii++){
	      forwardArray[i_anc][ii] = 0;
	      backwardArray[i_anc][ii] = 0;
	      for(int iii = 0; iii< nr_founders; iii++){
		forwardArray[i_anc][ii]  += forward_hap_array[ii][iii];
		backwardArray[i_anc][ii]  += backward_hap_array[ii][iii];
	      }
	      forwardArray[i_anc][ii] = log(forwardArray[i_anc][ii]);
	      backwardArray[i_anc][ii] = log(backwardArray[i_anc][ii]);
	      fprintf(stderr,"%d %d %f\n",i_anc,ii,forwardArray[i_anc][ii]);
	    }
	  }
	  
	  int best_i = 0;
	  double best_prob = -999999999;
	  for(int ii = 10; ii < 100; ii+=10){
	    if( forwardArray[0][ii]+backwardArray[1][ii] > best_prob){
	      best_i = ii;
	      best_prob = forwardArray[0][ii]+backwardArray[1][ii];
	    }
	  }
	  /*new breakpoint location*/
	  fprintf(stderr,"new location %d global (%d), prob (%f)\n", best_i,i_snp-50+best_i,best_prob);
	  
	  if(best_i < 50){
	    for(int ii = i_snp-50+best_i;ii<i_snp;ii++){
	      admx_samples_inf_anc[i_hap][ii] = anc[1]+'0';
	    }
	  }else{
	    for(int ii = i_snp; ii< i_snp-50+best_i; ii++){
	      admx_samples_inf_anc[i_hap][ii] = anc[0]+'0';
	    }
	  }
	  
	  i_snp+=50;
	}else {
	  i_snp++;
	}
      }
    }
  }
  char filename[100];
  sprintf(filename,"%s.phase",argv[8]);
  FILE *foutphase = fopen(filename,"w");
  sprintf(filename,"%s.anc",argv[8]);
  FILE *foutanc = fopen(filename,"w");
  /* print to output*/
  for(int i_hap = 0; i_hap < 4; i_hap++){
    for(int i_snp = 0 ; i_snp < nr_total_snps; i_snp++){
      fprintf(foutphase,"%c",admx_samples_inf_phase[i_hap][i_snp]);
      fprintf(foutanc,"%c",admx_samples_inf_anc[i_hap][i_snp]);
    }
    fprintf(foutphase,"\n");
    fprintf(foutanc,"\n");
  }
  fclose(foutphase);
  fclose(foutanc);
}






char* getgeno( char* hap1, char* hap2, int nr_snp){
	char *gen =  (char*)calloc (nr_total_snps+1,sizeof(char));
	for(int i =0; i < nr_snp; i++){
		gen[i] = (char)( (hap1[i] -48)+ (hap2[i]-48) +48);
	}
}


void read_true_sol(const char* filename){
	
	FILE* fin = fopen(filename,"r");
	char myhap[MAX_SNPS];
	int i =0;
	
	if(!fin){
		fprintf(stderr, "Could not open file: %s !!\n",filename);
		exit(1);
	}
	
	i = 0;
	while (fgets(myhap,MAX_SNPS,fin)) {
		i++;
	}
	assert (i==nr_admx_samples);
	
	admx_samples_anc = (char**)calloc(nr_admx_samples , sizeof(char*));	
	for(int i = 0; i < nr_admx_samples; i++)
	{
		admx_samples_anc[i] = (char*)calloc((nr_total_snps+1) , sizeof(char));
	}
	rewind(fin);	
	i=0;
	while (fgets(admx_samples_anc[i],MAX_SNPS,fin)) {
		i++;
	}
	fclose(fin);
	
}



void read_pos_file(const char* filename){
  FILE* fin = fopen(filename,"r");
  char myhap [MAX_SNPS];
  
  int i = 0;
  
  if(!fin){
    fprintf(stderr, "Could not open file: %s !!\n",filename);
    exit(1);
  }
  
  i = 0;
  while (fgets(myhap,MAX_SNPS,fin)) {
    i++;
  }
  
  nr_total_snps = i;
  fprintf(stderr,"reading %d snps from pos file\n",nr_total_snps);
  
  snps_pos = (int*)calloc(nr_total_snps, sizeof(int));
  
  rewind(fin);	
  i=0;
  while (fgets(myhap,MAX_SNPS,fin)) {
    snps_pos[i] = atol(myhap);
    i++;
  }
  fclose(fin);
}
void read_map_file(const char* filename){
	FILE* fin = fopen(filename,"r");
	char myhap [MAX_SNPS];
	
	int i = 0;
	
	if(!fin){
		fprintf(stderr, "Could not open file: %s !!\n",filename);
		exit(1);
	}
	
	i = 0;
	while (fgets(myhap,MAX_SNPS,fin)) {
		i++;
	}
	
	nr_total_snps = i;
	fprintf(stderr,"reading %d snps from map file\n",nr_total_snps);
	
	snps_map = (double*)calloc(nr_total_snps, sizeof(double));
	
	rewind(fin);	
	i=0;
	while (fgets(myhap,MAX_SNPS,fin)) {
		snps_map[i] = atof(myhap);
		i++;
	}
	fclose(fin);
}


void read_samples(const char* filename){
  FILE* fin = fopen(filename,"r");
  char m[MAX_SNPS];
  char f[MAX_SNPS];
  char c[MAX_SNPS];
  int i = 0;
  
  if(!fin){
    fprintf(stderr, "Could not open file: %s !!\n",filename);
    exit(1);
  }

  fgets(m,MAX_SNPS,fin);
  fgets(f,MAX_SNPS,fin);
  fgets(c,MAX_SNPS,fin);
  //create haplotype panels;
  //nr_total_snps = strlen(m);

  nr_admx_samples = 4;
  admx_samples = (char**)calloc(nr_admx_samples+1, sizeof(char*));	
  admx_samples_inf_anc = (char**)calloc(nr_admx_samples+1 , sizeof(char*));
  admx_samples_inf_phase = (char**)calloc(nr_admx_samples+1 , sizeof(char*));
  for(int i = 0; i < nr_admx_samples+1; i++)
    {
      admx_samples[i] = (char*)calloc((nr_total_snps+1) , sizeof(char));
      admx_samples_inf_anc[i] = (char*)calloc((nr_total_snps+1) , sizeof(char));
      admx_samples_inf_phase[i] = (char*)calloc(nr_total_snps+1 , sizeof(char));
    }
  for(int i = 0; i < nr_total_snps; i++){
    admx_samples[0][i] = admx_samples[1][i] = admx_samples[2][i] = admx_samples[3][i] = '?'; 

    if(m[i] == '0'){
      admx_samples[0][i] = admx_samples[1][i] = '0'; 
    }
    if(m[i] == '2'){
      admx_samples[0][i] = admx_samples[1][i] = '1'; 
    }
    if(f[i] == '0'){
      admx_samples[2][i] = admx_samples[3][i] = '0'; 
    }
    if(f[i] == '2'){
      admx_samples[2][i] = admx_samples[3][i] = '1'; 
    }
    if(c[i] == '0'){
      admx_samples[0][i] = admx_samples[2][i] = '0'; 
    }
    if(c[i] == '2'){
      admx_samples[0][i] = admx_samples[2][i] = '1'; 
    }
    for(int ii=0; ii< 3; ii++){
      if(m[i] == '1'){
	if(admx_samples[0][i]== '0'){
	  admx_samples[1][i] = '1';
	}
	if(admx_samples[0][i]== '1'){
	  admx_samples[1][i] = '0';
	}
      }
      if(f[i] == '1'){
	if(admx_samples[2][i]== '0'){
	  admx_samples[3][i] = '1';
	}
	if(admx_samples[2][i]== '1'){
	  admx_samples[3][i] = '0';
	}
      }
      if(c[i] == '1'){
	if(admx_samples[0][i] == '0'){
	  admx_samples[2][i] = '1';
	}
	if(admx_samples[0][i] == '1'){
	  admx_samples[2][i] = '0';
	}
	if(admx_samples[2][i] == '0'){
	  admx_samples[0][i] = '1';
	}
	if(admx_samples[2][i] == '1'){
	  admx_samples[0][i] = '0';
	}
      }
    }
    //fprintf(stderr,"%d (%c %c %c) - > (%c %c %c %c)",i,m[i],f[i],c[i],admx_samples[0][i],admx_samples[1][i],admx_samples[2][i],admx_samples[3][i]); getchar();
  }
}


void read_training_haps(int ancestry, const char* filename){
	
	FILE* fin = fopen(filename,"r");
	char myhap [MAX_SNPS];
	
	int i = 0;
	
	if(!fin){
		fprintf(stderr, "Could not open file: %s !!\n",filename);
		exit(1);
	}
	
	anc_map[ancestry].nr_hap = 0;
	while (fgets(myhap,MAX_SNPS,fin)) {
		anc_map[ancestry].nr_hap++;
	}
	
	anc_map[ancestry].hap_array = (char**)calloc(anc_map[ancestry].nr_hap , sizeof(char*));	
	for(int i = 0; i < anc_map[ancestry].nr_hap; i++)
	{
		anc_map[ancestry].hap_array[i] = (char*)calloc((nr_total_snps+1) , sizeof(char));
	}
	rewind(fin);	
	i=0;
	while (fgets(anc_map[ancestry].hap_array[i],MAX_SNPS,fin)) {
		i++;
	}
	fclose(fin);
}

/* updates hmm and returns true if anything is modifed, false otherwise*/
bool update_hmm(HAP_HMM *hmm,double* new_pi, STATE** new_states,double*** new_trans)
{
	/* make every prob higher than MIN_PROB*/
	/*1 pi*/
	for(int i =0; i< nr_founders; i++)
	{
		if( new_pi[i] < MIN_PROB)
			new_pi[i] = MIN_PROB;
	}
	/*2 transitions*/			
	for(int i =0; i < nr_snp-1; i++)
	{
		for(int j =0; j < nr_founders; j++)
		{
			
			for(int k =0; k < nr_founders; k++)
			{
				if( new_trans[i][j][k] < MIN_PROB)
					new_trans[i][j][k] = MIN_PROB;
				
			}
		}
	}
	/*3 emissions*/
	for(int i =0; i < nr_snp; i++)
	{
		for(int j =0; j < nr_founders; j++)
		{
			
			if(new_states[i][j].prob_em[0] < MIN_PROB)
				new_states[i][j].prob_em[0]  = MIN_PROB;
			if(new_states[i][j].prob_em[1] < MIN_PROB)
				new_states[i][j].prob_em[1]  = MIN_PROB;
			
		}
	}
	/* renormalize*/
	/*1 normalize pi*/
	double total=0.0;
	for(int i =0; i< nr_founders; i++)
	{
		total += new_pi[i];
	}
	for(int i =0; i< nr_founders; i++)
	{
		new_pi[i] /= total;
	}
	/*2 normalize transitions*/			
	
	for(int i =0; i < nr_snp-1; i++)
	{
		for(int j =0; j < nr_founders; j++)
		{
			total = 0.0;
			for(int k =0; k < nr_founders; k++)
			{
				total += new_trans[i][j][k];
				
			}
			for(int k =0; k < nr_founders; k++)
			{
				new_trans[i][j][k] /= total;
			}
		}
	}
	/*3 normalize emissions*/
	for(int i =0; i < nr_snp; i++)
	{
		for(int j =0; j < nr_founders; j++)
		{
			total = new_states[i][j].prob_em[0] + new_states[i][j].prob_em[1];
			new_states[i][j].prob_em[0] /= total;
			new_states[i][j].prob_em[1] /= total;				
		}
	}
	/* update HMM*/
	double is_modif = 0.0;
	double max_modif = 0.0;
	/*1 pi*/
	for(int i =0; i< nr_founders; i++)
	{
		if(fabs((*hmm).pi[i] - new_pi[i]) > max_modif)
			max_modif = fabs((*hmm).pi[i] - new_pi[i]);
		is_modif +=fabs((*hmm).pi[i] - new_pi[i]);		
		(*hmm).pi[i] = new_pi[i];
	}
	/*2 transitions*/			
	for(int i =0; i < nr_snp-1; i++)
	{
		for(int j =0; j < nr_founders; j++)
		{
			for(int k =0; k < nr_founders; k++)
			{
				if(fabs((*hmm).trans[i][j][k] - new_trans[i][j][k])> max_modif)
					max_modif = fabs((*hmm).trans[i][j][k] - new_trans[i][j][k]);
				is_modif += fabs((*hmm).trans[i][j][k] - new_trans[i][j][k]);
				(*hmm).trans[i][j][k] = new_trans[i][j][k] ;
				
			}
		}
	}
	/*3 emissions*/
	for(int i =0; i < nr_snp; i++)
	{
		for(int j =0; j < nr_founders; j++)
		{
			if(fabs((*hmm).states[i][j].prob_em[0] - new_states[i][j].prob_em[0])> max_modif)
				max_modif = fabs((*hmm).states[i][j].prob_em[0] - new_states[i][j].prob_em[0]);
			
			is_modif +=fabs((*hmm).states[i][j].prob_em[0] - new_states[i][j].prob_em[0]);
			(*hmm).states[i][j].prob_em[0] = new_states[i][j].prob_em[0];
			
			if(fabs((*hmm).states[i][j].prob_em[1] - new_states[i][j].prob_em[1])> max_modif)
				max_modif = fabs((*hmm).states[i][j].prob_em[1] - new_states[i][j].prob_em[1]);
			
			is_modif +=fabs((*hmm).states[i][j].prob_em[1] - new_states[i][j].prob_em[1]);
			(*hmm).states[i][j].prob_em[1] = new_states[i][j].prob_em[1];
		}
	}
	int nr_var;
	nr_var = nr_founders + (nr_snp-1)*nr_founders*nr_founders + 2*nr_founders*nr_snp;
	if(is_modif/nr_var > MIN_NR ) return true;
	return false;	
}
void alloc_mem_other()
{
	forward_hap_array = (double**)malloc(MAXBLOCK*sizeof(double*));
	for(int i =0; i < MAXBLOCK; i++)
	{
		forward_hap_array[i] = (double*)malloc(nr_founders*sizeof(double));
	}	
	backward_hap_array = (double**)malloc(MAXBLOCK*sizeof(double*));
	for(int i =0; i < MAXBLOCK; i++)
	{
		backward_hap_array[i] = (double*)malloc(nr_founders*sizeof(double));
	}
	
	//Uno
	array_uno_precomp1 = (double**)malloc(nr_founders*sizeof(double*));
	for(int i =0; i< nr_founders; i++)
	{
		array_uno_precomp1[i] = (double*)malloc(nr_founders*sizeof(double));
	}
	
	totalProbUnoMatrix_forward = (double***)malloc(MAXBLOCK*sizeof(double**)); 
	totalProbUnoMatrix_backward = (double***)malloc(MAXBLOCK*sizeof(double**));
	totalProbUnoMatrix_unk = (double**)malloc(MAXBLOCK*sizeof(double*));
	
	for(int i =0; i < MAXBLOCK; i++)
	{
		totalProbUnoMatrix_forward[i] = (double**)malloc(nr_founders*sizeof(double*)); //Genotype Total
		totalProbUnoMatrix_backward[i] = (double**)malloc(nr_founders*sizeof(double*)); //Genotype Total
		totalProbUnoMatrix_unk[i] = (double*)malloc(nr_founders*sizeof(double)); //Genotype Total
		for(int j =0; j < nr_founders; j++)
		{
			totalProbUnoMatrix_forward[i][j] = (double*)malloc(nr_founders*sizeof(double)); //Genotype Total
			totalProbUnoMatrix_backward[i][j] = (double*)malloc(nr_founders*sizeof(double)); //Genotype Total
		}
	}
	pi_temp = (double*) malloc (nr_founders*sizeof(double));
	trans_temp = (double***)malloc((MAXBLOCK-1)*sizeof(double**));
	for(int i =0; i < MAXBLOCK-1; i++)
	{
		trans_temp[i] = (double**)malloc(nr_founders*sizeof(double*));
		for(int j =0; j < nr_founders; j++)
		{
			trans_temp[i][j] = (double*)malloc(nr_founders*sizeof(double));
		}
	}
	
	states_temp = (STATE**)malloc(MAXBLOCK*sizeof(STATE*));
	for(int i =0; i < MAXBLOCK; i++)
	{
		states_temp[i] = (STATE*)malloc(nr_founders*sizeof(STATE));
	}
	
		
	
}
void alloc_mem_hmm(HAP_HMM *hmm)
{	
	hmm->states = (STATE**)malloc(MAXBLOCK*sizeof(STATE*));
	for(int i =0; i < MAXBLOCK; i++)
	{
		hmm->states[i] = (STATE*)malloc(nr_founders*sizeof(STATE));
	}
	hmm->pi = (double*)malloc(nr_founders*sizeof(double));
	hmm->pi_end = (double*)malloc(nr_founders*sizeof(double));
	hmm->trans = (double***)malloc((MAXBLOCK-1)*sizeof(double**));
	for(int i =0; i < MAXBLOCK-1; i++)
	{
		hmm->trans[i] = (double**)malloc(nr_founders*sizeof(double*));
		for(int j =0; j < nr_founders; j++)
		{
			hmm->trans[i][j] = (double*)malloc(nr_founders*sizeof(double));
		}
	}
}
void free_mem_other()
{
	
	for(int i =0; i < MAXBLOCK; i++)
		free(forward_hap_array[i]);
	free(forward_hap_array);
	
	for(int i =0; i < MAXBLOCK; i++)
		free(backward_hap_array[i]);
	free(backward_hap_array);
	//Uno
	for(int i =0; i< nr_founders; i++)
		free(array_uno_precomp1[i]);
	free(array_uno_precomp1);
	
	for(int i =0; i < MAXBLOCK; i++)
	{
		for(int j =0; j < nr_founders; j++)
		{
			free(totalProbUnoMatrix_forward[i][j]);
			free(totalProbUnoMatrix_backward[i][j]);
		}
		free(totalProbUnoMatrix_forward[i]);
		free(totalProbUnoMatrix_backward[i]);
		free(totalProbUnoMatrix_unk[i]);
	}
	free(totalProbUnoMatrix_forward);
	free(totalProbUnoMatrix_backward);
	free(totalProbUnoMatrix_unk);
}
void free_mem_hmm(HAP_HMM *hmm)
{
	for(int i =0; i < MAXBLOCK; i++)
		free(hmm->states[i]);
	free(hmm->states);
	free(hmm->pi);
	free(hmm->pi_end);
	for(int i =0; i < MAXBLOCK-1; i++)
	{
		for(int j =0; j < nr_founders; j++)
			free(hmm->trans[i][j]);
		free(hmm->trans[i]);
	}
	free(hmm->trans);
}
int compute_allele_occurrences(int snp,char alelle,int h, int start_pos)
{
	int freq = 0;
	for(int j = 0; j < anc_map[h].nr_hap; j++)
	{
		if(anc_map[h].hap_array[j][start_pos+snp] == alelle)
		{
			freq++;
		}
	}
	return freq;
}

void def_val_hmm(HAP_HMM *hmm, int h, int start_pos)
{
	double noise_dif = 0.0;
	for(int j = 0; j< nr_founders; j++)
	{
		hmm->pi[j] =  1/(double)nr_founders;
#if ADD_NOISE
		double p = pow(EPS, (unif_rand((long)(200*noise)) - 100.0*noise)/(double)100);
		noise_dif += fabs(hmm->pi[j] - hmm->pi[j] * p);
		hmm->pi[j] *= p;
#endif			
	}
	for(int j = 0; j< nr_founders; j++)
		hmm->pi_end[j] =  1/(double)nr_founders;
	
	
	
	/*set transition to 1/nr_founders for every trans*/	
	for(int i = 0; i< nr_snp-1; i++)
	{
		for(int j = 0; j< nr_founders; j++)
		{
			for(int k = 0; k< nr_founders; k++)
			{
				if(k==j)
				{
					hmm->trans[i][j][k] = founder_bias;
				}else
				{
					hmm->trans[i][j][k] = (1-founder_bias)/(double)(nr_founders-1);
				}
#if ADD_NOISE			
				double p = pow(EPS, (unif_rand((long)(200*noise)) - 100.0*noise)/(double)100);
				noise_dif += fabs(hmm->trans[i][j][k]  - hmm->trans[i][j][k]  * p);
				hmm->trans[i][j][k]  *= p;
#endif			
			}
		}
	}	
	/* set emission prob according to alele freq*/
	for(int i = 0; i< nr_snp; i++)
	{
		int num_0 = compute_allele_occurrences(i,'0',h,start_pos) + 1;
		int num_1 = compute_allele_occurrences(i,'1',h,start_pos) + 1;
#ifndef SAME_EM_PROB
		int tot = num_0 + num_1;
		for(int j = 0; j< nr_founders; j++)
		{
			if(unif_rand(tot) < num_0)
			{
				hmm->states[i][j].major = 0;
				hmm->states[i][j].prob_em[0] = MAJOR_PROB;
				hmm->states[i][j].prob_em[1] = 1 - MAJOR_PROB;					
			}else
			{
				hmm->states[i][j].major = 1;
				hmm->states[i][j].prob_em[0] = 1 - MAJOR_PROB;					
				hmm->states[i][j].prob_em[1] = MAJOR_PROB;
			}
		}
#endif
#ifdef SAME_EM_PROB
		if(num_0 == num_1)
		{
			for(int j = 0; j< nr_founders; j++)
				hmm->states[i][j].prob_em[0] = hmm->states[i][j].prob_em[1] = 0.5;
		}
		else
		{
			for(int j = 0; j< nr_founders; j++)
			{
				hmm->states[i][j].prob_em[0] = num_0/double(num_0+num_1) ;
				hmm->states[i][j].prob_em[1] = num_1/double(num_0+num_1) ;
			}	
		}
#endif		
#if ADD_NOISE	
		for(int j = 0; j< nr_founders; j++)
		{
			double p = pow(EPS, (unif_rand((long)(200*noise)) - 100.0*noise)/(double)100);
			noise_dif += fabs(hmm->states[i][j].prob_em[0]  - hmm->states[i][j].prob_em[0]  * p);
			hmm->states[i][j].prob_em[0] *= p;
			
			p = pow(EPS, (unif_rand((long)(200*noise)) - 100.0*noise)/(double)100);
			noise_dif += fabs(hmm->states[i][j].prob_em[1]  - hmm->states[i][j].prob_em[1]  * p);
			hmm->states[i][j].prob_em[1] *= p;
			
		}		
		
#endif			
	}
#if ADD_NOISE		
	/*if noise was added then renormalize*/	
	/*1 normalize pi*/
	double total=0.0;
	for(int i =0; i< nr_founders; i++)
	{
		total += hmm->pi[i];
	}
	for(int i =0; i< nr_founders; i++)
	{
		hmm->pi[i] /= total;
		
	}
	/*2 normalize transitions*/			
	
	for(int i =0; i < nr_snp-1; i++)
	{
		for(int j =0; j < nr_founders; j++)
		{
			total = 0.0;
			for(int k =0; k < nr_founders; k++)
			{
				total += hmm->trans[i][j][k];
				
			}
			for(int k =0; k < nr_founders; k++)
			{
				hmm->trans[i][j][k] /= total;
			}
		}
	}
	/*3 normalize emissions*/
	for(int i =0; i < nr_snp; i++)
	{
		for(int j =0; j < nr_founders; j++)
		{
			total = hmm->states[i][j].prob_em[0] + hmm->states[i][j].prob_em[1];
			hmm->states[i][j].prob_em[0] /= total;
			hmm->states[i][j].prob_em[1] /= total;				
		}
	}
#endif						
}
double forward_hap(int a, int j , char* hap ,HAP_HMM hmm, int start_pos)
{
	double ret;
	if(j == 0)
	{
		ret = hmm.pi[a];
		if((hap[start_pos+j]== '0')) 
			ret *= hmm.states[j][a].prob_em[0]; 
		else if((hap[start_pos+j]== '1'))	
			ret *= hmm.states[j][a].prob_em[1]; 
		else if((hap[start_pos+j]== '?')||(hap[start_pos+j]=='3'))
			ret *= 1.0;
		else
		{
			fprintf(stderr,"1ERROR forward_hap hap[%d]=%c not found\n\n", start_pos+j,hap[start_pos+j] );
			exit(0);
		}
		return ret;
		
	}
	else
	{
		double tot = 0.0;
		for(int prev_a = 0 ; prev_a < nr_founders; prev_a++)
		{
			double prob = hmm.trans[j-1][prev_a][a] * forward_hap_array[j-1][prev_a];
			if((hap[start_pos+j]== '0'))
				prob *= hmm.states[j][a].prob_em[0]; 
			else if((hap[start_pos+j]== '1'))
				prob *= hmm.states[j][a].prob_em[1]; 
			else if((hap[start_pos+j]== '?')||(hap[start_pos+j]=='3')) 
				prob *=1.0;
			else
			{
				fprintf(stderr,"2ERROR forward_hap hap[%d]=%c not found\n\n", start_pos+j,hap[start_pos+j] );
				exit(0);
			}
			tot += prob;			
		}
		return tot;
	}
}

double backward_hap(int a, int j , char* hap ,HAP_HMM hmm, int start_pos)
{
	
	if(j == nr_snp-1 ) return 1;
	else
	{
		double tot = 0.0;
		for(int next_a = 0 ; next_a < nr_founders; next_a++)
		{
			double prob = hmm.trans[j][a][next_a] * backward_hap_array[j+1][next_a];
			if((hap[start_pos+j+1]== '0')) 
				prob *= hmm.states[j+1][next_a].prob_em[0]; 
			else if((hap[start_pos+j+1]== '1')) 
				prob *= hmm.states[j+1][next_a].prob_em[1]; 
			else if((hap[start_pos+j+1]== '?')||(hap[start_pos+j+1]=='3')) 
				prob *= 1.0;
			else
			{
				fprintf(stderr,"ERROR backward_hap hap[%d]=%c not found\n\n", start_pos+j,hap[start_pos+j] );
				exit(0);
			}
			tot += prob;			
		}
		return tot;
	}
}
double gethapprob(char* hap ,HAP_HMM hmm, int start_pos){
	for(int i =0 ;i<nr_snp;i++)
	{
		for(int a = 0 ; a < nr_founders; a++)
		{
			forward_hap_array[i][a] = forward_hap( a, i , hap , hmm,start_pos);
		}
	}
	double ret = 0;
	for(int a = 0 ; a < nr_founders; a++){
		ret += forward_hap_array[nr_snp-1][a];
	}
	return (ret);
	
}

void fill_forward_backward_hap_arrays(char* hap ,HAP_HMM hmm, int start_pos)
{
	for(int i =0 ;i<nr_snp;i++)
	{
		for(int a = 0 ; a < nr_founders; a++)
		{
			forward_hap_array[i][a] = forward_hap( a, i , hap , hmm,start_pos);
		}
	}
	for(int i =nr_snp-1 ;i>=0;i--)
	{
		for(int a = 0 ; a < nr_founders; a++)
		{
			backward_hap_array[i][a] = backward_hap( a, i , hap , hmm,start_pos);
		}
	}
}
double getMinInitDenom(int h)//set numerator=denominator so transitions are 1 probability for initialization
{
	if(nr_founders > 1) return (nr_founders-1.0);
	else return TRANSITION_PSEUDOCOUNT*anc_map[h].nr_hap;
}




void trainHMM(int h,int start_pos)
{
	
		
	//fprintf(stderr,"Starting training HMM at anc %d\n",h);
	
	def_val_hmm(&anc_map[h].hmm,h,start_pos);
	
	

	bool cont = true;
	int step =0;
	double prev_prob=- 1/MIN_NR;
	while (cont)
	{	
		/*initialize temp structures*/
		for(int i =0; i < nr_founders; i++)
			pi_temp[i] = TRANSITION_PSEUDOCOUNT*anc_map[h].nr_hap/getMinInitDenom(h); /* don't let any prob get 0*/
		for(int i =0; i < nr_snp-1; i++)
			for(int j =0; j < nr_founders; j++)
				for(int k =0; k < nr_founders; k++)
					trans_temp[i][j][k] = TRANSITION_PSEUDOCOUNT*anc_map[h].nr_hap/getMinInitDenom(h); /* don't let any prob get 0*/
		for(int i =0; i < nr_snp; i++)
			for(int j =0; j < nr_founders; j++)
			{
		 		states_temp[i][j].prob_em[0] = EMISSION_PSEUDOCOUNT*anc_map[h].nr_hap/nr_founders; /* don't let any prob get 0*/
	 			states_temp[i][j].prob_em[1] = EMISSION_PSEUDOCOUNT*anc_map[h].nr_hap/nr_founders; /* don't let any prob get 0*/
			}
		double tot_prob = 0.0;
		double p_paths= 0.0;
		double p_haps= 0.0;
		double prob_of_hap;
		// use baum-welch training
		tot_prob = 0.0;
		for(int index_hap=0;  index_hap < anc_map[h].nr_hap;index_hap++)
		{
 			//compute forward_backward
		  fill_forward_backward_hap_arrays(anc_map[h].hap_array[index_hap] , anc_map[h].hmm,start_pos);
			//compute prob of hap
			prob_of_hap=0.0;
			for(int i = 0 ; i < nr_founders; i++)
				prob_of_hap += forward_hap_array[nr_snp-1][i];		
			tot_prob += log(prob_of_hap);
			for(int t = 0 ; t < nr_snp-1; t++)
			{
		 		//update trans_temp
		 		for(int i = 0 ; i < nr_founders; i++)
				{
					for(int j = 0 ; j < nr_founders; j++)
					{
				   		if(forward_hap_array[t][i] *
						   anc_map[h].hmm.trans[t][i][j] *
						   anc_map[h].hmm.states[t+1][j].prob_em[ anc_map[h].hap_array[index_hap][start_pos+ t+1]-48 ]*
						   backward_hap_array[t+1][j] / prob_of_hap < MIN_NR)
						{
							fprintf(stderr,"ERROWRRRRRRRR!!! %e %e\n",MIN_NR, forward_hap_array[t][i] *
									anc_map[h].hmm.trans[t][i][j] *
									anc_map[h].hmm.states[t+1][j].prob_em[ anc_map[h].hap_array[index_hap][start_pos+t+1]-48 ]*
									backward_hap_array[t+1][j] / prob_of_hap);
						}					
						trans_temp[t][i][j] += forward_hap_array[t][i] *
						anc_map[h].hmm.trans[t][i][j] *
						anc_map[h].hmm.states[t+1][j].prob_em[ anc_map[h].hap_array[index_hap][start_pos+t+1]-48 ]*
						backward_hap_array[t+1][j] / prob_of_hap ;
					}
				}
	 		}
			for(int t = 0 ; t < nr_snp; t++)		
				for(int i =0; i < nr_founders; i++)
					states_temp[t][i].prob_em[ anc_map[h].hap_array[index_hap][start_pos+t]-48 ] += forward_hap_array[t][i] * backward_hap_array[t][i] / prob_of_hap;
			for(int i =0; i < nr_founders; i++)
				pi_temp[i] += anc_map[h].hmm.states[0][i].prob_em[ anc_map[h].hap_array[index_hap][start_pos+0]-48 ]*backward_hap_array[0][i] / prob_of_hap;
		}
		//end training
		/*make all temp structures probabilities and update hmm*/
		/*1 normalize pi*/
		double total=0.0;
		for(int i =0; i< nr_founders; i++)
			total += pi_temp[i];
		for(int i =0; i< nr_founders; i++)
			pi_temp[i] /= total;
		/*2 normalize transitions*/			
		for(int i =0; i < nr_snp-1; i++)
			for(int j =0; j < nr_founders; j++)
			{
				total = 0.0;
				for(int k =0; k < nr_founders; k++)
					total += trans_temp[i][j][k];
				for(int k =0; k < nr_founders; k++)
					trans_temp[i][j][k] /= total;
			}
		/*3 normalize emissions*/
		for(int i =0; i < nr_snp; i++)
			for(int j =0; j < nr_founders; j++)
			{
				total = states_temp[i][j].prob_em[0] + states_temp[i][j].prob_em[1];
				states_temp[i][j].prob_em[0] /= total;
				states_temp[i][j].prob_em[1] /= total;				
			}
		cont = update_hmm(&anc_map[h].hmm,pi_temp,states_temp,trans_temp);	
		step++;
		if(step == max_steps)
			cont = false;
		//Baum Welch stopping condition
		//fprintf(stderr,"Percentage increase %f(%f)\n", 1 - tot_prob/prev_prob);
		//if(tot_prob - prev_prob < 3)
		if( 1 - tot_prob/prev_prob < 0.01)
			cont = false;
		prev_prob = tot_prob;	
	}//end training
	
	//fprintf(stderr,"---- %d steps\n",step);
}


int convAlleleCharToInteger(char val)
{
	if(val == '0') return 0;
	else if (val == '1') return 1;
	else if (val == '2') return 2;
	else if (val == '3') return 3;
	else return 3;//unknown
}
int convAlleleIntegerToChar(int val)
{
	if(val == 0) return '0';
	else if (val == 1) return '1';
	else if (val == 2) return '2';
	else if (val == 3) return '3';
	else return '3';//unknown
}
double calcEmissionProbMax(int g1, int g2, char ch,int j, HAP_HMM hmmOne, HAP_HMM hmmTwo)
//given the genotype at location SNP j, find the hmm probability
//that it will emit a compatible set of characters, at the states of the 2 paths g1,g2
{
	double maxVal = 0;
	double currentVal = 0;
	double g10 = hmmOne.states[j][g1].prob_em[0];
	double g11 = hmmOne.states[j][g1].prob_em[1];
	double g20;
	double g21;
	g20 = hmmTwo.states[j][g2].prob_em[0];
	g21 = hmmTwo.states[j][g2].prob_em[1];
	int chi = convAlleleCharToInteger(ch);
	
	currentVal = g10*g20;
	if (((chi==0)||(chi==3))&&(currentVal>maxVal))
		maxVal = currentVal;
	currentVal = g10*g21;
	if (((chi==1)||(chi==3))&&(currentVal>maxVal))
		maxVal = currentVal;
	currentVal = g11*g20;
	if (((chi==1)||(chi==3))&&(currentVal>maxVal))
		maxVal = currentVal;
	currentVal = g11*g21;
	if (((chi==2)||(chi==3))&&(currentVal>maxVal))
		maxVal = currentVal;
	return maxVal;
}



double totalProbUnoAt_forward_fast(int a, int b, int j, int chi, HAP_HMM hmmOne, HAP_HMM hmmTwo)
{
	double ret = 0;
	double ep = calcEmissionProbMax(a,b,convAlleleIntegerToChar(chi),j,hmmOne,hmmTwo);
	if(j == 0) 
	{
		ret = hmmOne.pi[a]*hmmTwo.pi[b]*ep;		
		return ret;
	}
	else
	{
		for(int prev_b = 0 ; prev_b < nr_founders; prev_b++)
			ret+= array_uno_precomp1[a][prev_b]*hmmTwo.trans[j-1][prev_b][b];
		return ret*ep;
	}
}

double preCompute_total_uno1_forward(int a, int bpr, int j , HAP_HMM hmm)
{
	double ret = 0;
	for(int prev_a = 0 ; prev_a < nr_founders; prev_a++)
		ret += totalProbUnoMatrix_forward[j-1][prev_a][bpr]*hmm.trans[j-1][prev_a][a];
	return ret;
}
double totalProbUno_forward_fast(char* ch, HAP_HMM hmmOne, HAP_HMM hmmTwo,int start_pos)
//Genotype Total...Given the input haplotype and data hidden markov model, compute the summed probability that this haplotype will be emitted
{
	double ret = 0.0;
	int chi;
	for (int j = 0; j < nr_snp;j++)
	{
		chi = convAlleleCharToInteger(ch[start_pos+j]);
		for (int a =0; a< nr_founders;a++)
		{
			if(j>0)
			{	
				for (int b =0; b< nr_founders;b++)
					array_uno_precomp1[a][b] = preCompute_total_uno1_forward(a,b,j,hmmOne);
			}
			for (int b =0; b< nr_founders;b++)
			{
			 	totalProbUnoMatrix_forward[j][a][b] = totalProbUnoAt_forward_fast(a,b,j,chi,hmmOne,hmmTwo);
				if (j == nr_snp-1) 
					ret+=totalProbUnoMatrix_forward[j][a][b]*hmmOne.pi_end[a]*hmmTwo.pi_end[b];
			}
		}
		//fprintf(stderr,"%d %e\n",j,totalProbUnoMatrix_forward[j][0][0] );
	}
	return ret;
}



double preCompute_total_uno1_backward(int a, int bpr, int j , HAP_HMM hmm)
{
	double ret = 0;
	for(int prev_a = 0 ; prev_a < nr_founders; prev_a++)
		ret += totalProbUnoMatrix_backward[j+1][prev_a][bpr]*hmm.trans[j][a][prev_a];
	return ret;
}
double totalProbUnoAt_backward_fast(int a, int b, int j, int chi, HAP_HMM hmmOne, HAP_HMM hmmTwo)
{
	double ret = 0;
	double ep = calcEmissionProbMax(a,b,convAlleleIntegerToChar(chi),j,hmmOne, hmmTwo);
	if(j == (nr_snp-1)) return hmmOne.pi_end[a]*hmmTwo.pi_end[b]*ep;		
	else
	{
		for(int prev_b = 0 ; prev_b < nr_founders; prev_b++)
			ret+= array_uno_precomp1[a][prev_b]*hmmTwo.trans[j][b][prev_b];
		return ret*ep;
	}
}

double totalProbUno_backward_fast(char* ch, HAP_HMM hmmOne, HAP_HMM hmmTwo,int start_pos)
//Genotype Total...Given the input haplotype and data hidden markov model, compute the summed probability that this haplotype set will be emitted
{
	double ret = 0;
	int chi;
	for (int j = (nr_snp-1); j >= 0;j--)
	{
		chi = convAlleleCharToInteger(ch[start_pos+j]);
		for (int a =0; a< nr_founders;a++)
		{
			if(j<(nr_snp-1))
			{	
				for (int b =0; b< nr_founders;b++)
					array_uno_precomp1[a][b] = preCompute_total_uno1_backward(a,b,j,hmmOne);
			}
			for (int b =0; b< nr_founders;b++)
			{
			 	totalProbUnoMatrix_backward[j][a][b] = totalProbUnoAt_backward_fast(a,b,j,chi,hmmOne, hmmTwo);
				if (j == 0) 
					ret+=totalProbUnoMatrix_backward[j][a][b]*hmmOne.pi[a]*hmmTwo.pi[b];
			}
			
		}
	}
	return ret;
}
