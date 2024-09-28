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
/* defs*/
#define MAX_SNPS 55000
#define MAJOR_EM_PROB 0.9
#define ADD_NOISE 1
#define LOGZERO -100000 

//#include <gsl/gsl_statistics_double.h>

double totalProbUno_forward_fast(char* ch, HAP_HMM hmmOne, HAP_HMM hmmTwo,int start_pos);
double totalProbUno_backward_fast(char* ch, HAP_HMM hmmOne, HAP_HMM hmmTwo,int start_pos);
/*variables*/
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
char **admx_samples_inf;
double ****admx_inf_array;


double eps = 0.1;
double founder_bias = .9;
int max_steps = 100;
int hmmCount = 5;
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

/*functions*/
void fill_forward_backward_hap_arrays(char* hap ,HAP_HMM hmm,int start_pos);
void read_training_haps(int ancestry, const char* filename);
void read_samples(const char* filename);
void read_true_sol(const char* filename);
void read_pos_file(const char* filename);
void read_map_file(const char* filename);
void trainHMM(int h,int start_pos);
void alloc_mem_other();
void get_tot_acc(int, FILE*);
int convAlleleCharToInteger(char val);
double mylog(double x){if(x<=0) return 0; return  log(x);};
void fix_ancestry(int right_snp);
void fix_ancestry_hap(int right_snp);
char* getgeno(char* hap1, char* hap2, int nr_snp);

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


int starting_snp = 0;
double *pi_temp;
double*** trans_temp;
STATE**	states_temp;



int main(int argc, char** argv){
	
	double alpha[5];
	int tot_c,tot_e;
	//	double anc_gen_forward[3][3][MAX_SNPS];
	//double anc_gen_backward[3][3][MAX_SNPS];
	
	
	double ****forwardUnoProb;
	double ****unoProb;
	int ****backPointerI;
	int ****backPointerJ;
	double bkpt_prob[MAX_SNPS];
	int best_s = -1;
	double best_l;
	int left[2],right[2],prev_win_r[2];
	int best_no_bkpt[2];
	int start_pos = 0;
	
	
	int lastwin =0;
	int winsize;
	int increment;
	
	nr_founders = 10;
	
	nr_anc = 5;
	hmmCount = nr_anc;
	recomb_rate =1e-8;
	num_gen = 15;
	
	init_rand(1);
	anc_map  = (ANC_HMM*)calloc(nr_anc, sizeof(ANC_HMM));
	
	if(argc != 11){
	  fprintf(stderr,"USAGE:: %s <win-size> <nr-founders> <SNP-file> <ANC1-haps> <ANC2-haps> <ANC3-haps>  <ANC4-haps> <ANC5-haps> <ADMX-geno> <out-file>\n", argv[0]);
	  exit(0);
	}
	
	//read_pos_file("/n/soco/xv/bioinfo/x/bogdan/lamp4/hapmap2/hapmix/HapmixRelease/jpt-chb-new/chr1.pos");
	nr_snp = atol(argv[1]);
	nr_founders = atol(argv[2]);
	read_pos_file(argv[3]);
	read_training_haps(0,argv[4]);
	read_training_haps(1,argv[5]);
	read_training_haps(2,argv[6]);
	read_training_haps(3,argv[7]);
	read_training_haps(4,argv[8]);
	read_samples(argv[9]);
	//read_true_sol("haps.anc");
	//read_map_file("/home/bp45/projects/localanc/latino-trio/latino/chr1.genmap");
	FILE *fout = fopen(argv[10],"w");
	
	window_padding = 50;
	increment = 50; 
	
	
	fprintf(stderr,"Read %d anc0 %d anc1 %d anc2 training haplotypes\n",anc_map[0].nr_hap,anc_map[1].nr_hap,anc_map[2].nr_hap);
	fprintf(stderr,"Read %d admixed samples\n",nr_admx_samples);
	fprintf(stderr,"Using window size %d with padding %d\n",nr_snp, window_padding);
	
	
	//set admx_inf_array
	for(int i =0; i < nr_admx_samples;i++){
		for(int j =0; j < nr_total_snps;j++){
			for(int ii =0; ii < nr_anc;ii++){
				for(int jj =0; jj < nr_anc; jj++){
					admx_inf_array[i][j][ii][jj] = 0;
				}
			}
		}
	}
	
	winsize = (int)(sqrt(eps)/(2*(num_gen-1)*recomb_rate*(1-alpha[0]*alpha[0])*(1-alpha[1]*alpha[1])));
	
	fprintf(stderr,"Using window size %d bp\n",winsize);
	
	forwardUnoProb = (double****)calloc(nr_total_snps/nr_snp,sizeof(double***));
	unoProb = (double****)calloc(nr_total_snps/nr_snp,sizeof(double***));
	backPointerI = (int****)calloc(nr_total_snps/nr_snp,sizeof(int***));
	backPointerJ = (int****)calloc(nr_total_snps/nr_snp,sizeof(int***));
	
	for(int i = 0; i < nr_total_snps/nr_snp; i++){
		forwardUnoProb[i] = (double***)calloc(nr_admx_samples,sizeof(double**));
		unoProb[i] = (double***)calloc(nr_admx_samples,sizeof(double**));
		backPointerI[i] = (int***)calloc(nr_admx_samples,sizeof(int**));
		backPointerJ[i] = (int***)calloc(nr_admx_samples,sizeof(int**));
		for(int j = 0; j < nr_admx_samples; j++){
			forwardUnoProb[i][j] = (double**)calloc(nr_anc,sizeof(double*));
			unoProb[i][j] = (double**)calloc(nr_anc,sizeof(double*));
			backPointerI[i][j] = (int**)calloc(nr_anc,sizeof(int*));
			backPointerJ[i][j] = (int**)calloc(nr_anc,sizeof(int*));
			for(int k = 0; k < nr_anc; k++){
				forwardUnoProb[i][j][k] = (double*)calloc(nr_anc,sizeof(double));
				unoProb[i][j][k] = (double*)calloc(nr_anc,sizeof(double));
				backPointerI[i][j][k] = (int*)calloc(nr_anc,sizeof(int));
				backPointerJ[i][j][k] = (int*)calloc(nr_anc,sizeof(int));
			}
		}

	}
	
	
	pi_temp = (double*) malloc (nr_founders*sizeof(double));
	trans_temp = (double***)malloc((nr_snp-1)*sizeof(double**));
	for(int i =0; i < nr_snp-1; i++)
	{
 		trans_temp[i] = (double**)malloc(nr_founders*sizeof(double*));
		for(int j =0; j < nr_founders; j++)
		{
			trans_temp[i][j] = (double*)malloc(nr_founders*sizeof(double));
		}
	}
	states_temp = (STATE**)malloc(nr_snp*sizeof(STATE*));
	for(int i =0; i < nr_snp; i++)
	{
 		states_temp[i] = (STATE*)malloc(nr_founders*sizeof(STATE));
	}
	
	for(int h = 0; h < hmmCount; h++){
		alloc_mem_hmm(&anc_map[h].hmm);
	}
	
	alloc_mem_other();
	
	for(int h = 0; h < hmmCount; h++){
		anc_map[h].forward_hap_array = (double**)malloc(nr_snp*sizeof(double*));
		for(int i =0; i < nr_snp; i++){
			anc_map[h].forward_hap_array[i] = (double*)malloc(nr_founders*sizeof(double));
		}	
		anc_map[h].backward_hap_array = (double**)malloc(nr_snp*sizeof(double*));
		for(int i =0; i < nr_snp; i++)    {
			anc_map[h].backward_hap_array[i] = (double*)malloc(nr_founders*sizeof(double));
		}
	}
	

	fprintf(stderr,"nr_total_snps %d\n",nr_total_snps);
	for(int idx_win = 0; idx_win < (nr_total_snps/nr_snp); idx_win++){
		start_pos = idx_win*nr_snp;
		fprintf(stderr,"computing at window %d [%d %d]\n",idx_win,idx_win*nr_snp,idx_win*nr_snp+nr_snp-1);
		for(int h = 0; h < hmmCount; h++){
			trainHMM(h,start_pos);
		}
		//fprintf(stderr,"training done\n");
		for(int igen = 0; igen < nr_admx_samples; igen++){
			for(int i =0; i < nr_anc; i++){
				for(int j = i ; j < nr_anc; j++){
				  unoProb[idx_win][igen][i][j]  = log(totalProbUno_forward_fast(admx_samples[igen], anc_map[i].hmm,anc_map[j].hmm,start_pos));
				  //fprintf(stderr,"%d %d (%d %d) %f\n",igen,idx_win,i,j,unoProb[idx_win][igen][i][j]);					
				}
			}
		}
		//fprintf(stderr,"training done\n");
	}
	
	
	for(int idx_win = 0; idx_win < (nr_total_snps/nr_snp); idx_win++){
	  if(idx_win == 0){
	    for(int igen = 0; igen < nr_admx_samples; igen++){
	      for(int i =0; i < nr_anc; i++){
		for(int j = i ; j < nr_anc; j++){
		  forwardUnoProb[idx_win][igen][i][j] = unoProb[idx_win][igen][i][j];
		  fprintf(stderr,"%d %d (-1 -1) ->(%d %d) %f\n",igen,idx_win,i,j,forwardUnoProb[idx_win][igen][i][j]);						
		}
	      }
	    }
	  }else{
	    for(int igen = 0; igen < nr_admx_samples; igen++){
	      
	      for(int i =0; i < nr_anc; i++){
		for(int j = i ; j < nr_anc; j++){
		  double maxp = -999999;
		  int max_i;
		  int max_j;
		  for(int prev_i =0; prev_i < nr_anc; prev_i++){
		    for(int prev_j = prev_i ; prev_j < nr_anc; prev_j++){
		      double trans = 1;	
		      if(prev_i != i ){
			trans = recomb_rate;//(snps_map[idx_win*nr_snp] - snps_map[(idx_win*nr_snp])/100;
		      }
		      if(prev_j != j ){
			trans  = recomb_rate;//(snps_map[idx_win*nr_snp] - snps_map[(idx_win*nr_snp])/100;
		      }
		      double prob  = log(trans) + (forwardUnoProb[idx_win-1][igen][prev_i][prev_j])
			+ (unoProb[idx_win][igen][i][j]);
		      if(prob > maxp){
			maxp = prob;
			max_i = prev_i;
			max_j = prev_j;
		      }
		      //fprintf(stderr," %d %d -> %d %d , %e %e\n",prev_i,prev_j,i,j, trans, unoProb[idx_win-1][igen][prev_i][prev_j]);
		    }
		  }
		  forwardUnoProb[idx_win][igen][i][j] = maxp;
		  backPointerI[idx_win][igen][i][j] = max_i;
		  backPointerJ[idx_win][igen][i][j] = max_j;
		  fprintf(stderr,"%d %d (%d %d) ->(%d %d) %f\n",igen,idx_win,max_i,max_j,i,j,forwardUnoProb[idx_win][igen][i][j]);						
		}
	      }
	    }
	  }
	}
	for(int igen = 0; igen < nr_admx_samples; igen++){
	  
	  double max = -99999;
	  int idx_win = (nr_total_snps/nr_snp)-1;
	  int cur_i = 0;
	  int cur_j = 0;
	  for(int i =0; i < nr_anc; i++){
	    for(int j = i ; j < nr_anc; j++){
	      if(forwardUnoProb[idx_win][igen][i][j] > max){
		max = forwardUnoProb[idx_win][igen][i][j];
		cur_i = i;
		cur_j = j;
	      }
	    }
		}
	  fprintf(stderr,"found max at %f (%d %d)\n",max,cur_i,cur_j);
	  while(idx_win >= 0){			
	    fprintf(stderr,"%d (%d %d)\n",idx_win,cur_i,cur_j);
	    
	    for(int i = idx_win*nr_snp; i<(idx_win+1)*nr_snp && i< nr_total_snps ;i++){
	      admx_inf_array[igen][i][cur_i][cur_j]++;
	    }
	    //fill in last bits of ancestry at the end of chromosome
	    if(nr_total_snps - (idx_win+1)*nr_snp <= nr_snp){
	      for(int i = (idx_win+1)*nr_snp; i< nr_total_snps ;i++){
		admx_inf_array[igen][i][cur_i][cur_j]++;
	      }
	    }
			
	    int tcur_i = backPointerI[idx_win][igen][cur_i][cur_j];
	    int tcur_j = backPointerJ[idx_win][igen][cur_i][cur_j];
	    //fprintf(stderr,"%d %d\n",tcur_i,tcur_j);
	    cur_i = tcur_i;
	    cur_j = tcur_j;
	    idx_win--;
	    
	  }
	}

	fix_ancestry_hap(nr_total_snps);
	/*
	  for(int igen = 0 ; igen< nr_admx_samples; igen++){		
	  for(int i =0; i < nr_total_snps; i++){
	  fprintf(fout," %d",admx_samples_inf[igen][i]);
	  }
	  fprintf(fout,"\n");
	  }
	*/
	
	int ambgPos=0;
	double forwardArray[1000];
	double backwardArray[1000];
	int anc[2];
	//nr_snp = 00;                                                                                                                                     


	for(int i_snp = nr_snp; i_snp < (nr_total_snps-nr_snp); i_snp+=nr_snp){
	  //fprintf(stderr,"%d %f %f %f\n",i_snp-1, getavganc(i_snp-1,0),getavganc(i_snp-1,1),getavganc(i_snp-1,2));                                       
	  //fprintf(stderr,"%d %f %f %f\n",i_snp, getavganc(i_snp,0),getavganc(i_snp,1),getavganc(i_snp,2));                                               
	  for(int h = 0; h < hmmCount; h++){
	    trainHMM(h, i_snp - (nr_snp/2));
	  }
	  for(int i_hap = 0; i_hap < nr_admx_samples; i_hap++){
	    if(admx_samples_inf[i_hap][i_snp] != admx_samples_inf[i_hap][i_snp-1]){

	      /*found a breakpoint -> readjusting it*/
	      int left[2];
	      int right[2];

	      left[0]= (int)(admx_samples_inf[i_hap][i_snp-1]/10);
	      left[1]= (int)(admx_samples_inf[i_hap][i_snp-1]%10);
	      right[0]= (int)(admx_samples_inf[i_hap][i_snp]/10);
	      right[1]= (int)(admx_samples_inf[i_hap][i_snp]%10);

	      //fprintf(stderr,"%d %d ( %2.0d[%d %d] -> %2.0d[%d,%d] )\n",i_hap,i_snp, admx_samples_inf[i_hap][i_snp-1],left[0], left[1], admx_samples_inf[i_hap][i_snp],right[0], right[1]);                     

	      totalProbUno_forward_fast(admx_samples[i_hap], anc_map[left[0]].hmm, anc_map[left[1]].hmm, i_snp-(nr_snp/2) );
	      totalProbUno_backward_fast(admx_samples[i_hap], anc_map[right[0]].hmm, anc_map[right[1]].hmm, i_snp-(nr_snp/2) );

	      for(int ii = 0; ii < nr_snp; ii++){
		forwardArray[ii] = 0;
		backwardArray[ii] = 0;
		for(int a = 0; a < nr_founders; a++){
		  for(int b = 0; b < nr_founders; b++){
		    forwardArray[ii]  += totalProbUnoMatrix_forward[ii][a][b];
		    backwardArray[ii] += totalProbUnoMatrix_backward[ii][a][b];
		  }
		}
		forwardArray[ii] = log(forwardArray[ii]);
		backwardArray[ii] = log(backwardArray[ii]);
		//fprintf(stderr,"%d %f %f\n",ii,forwardArray[ii],backwardArray[ii]);                                                                     
	      }
	      

	      int best_i = nr_snp/2;
	      double best_prob = MIN_LOG;
	      for(int ii = 10; ii < nr_snp-10; ii++){
		double rprob = 0;
		if(ii>0){
		  //rprob = log(recomb_rate*(snps_pos[i_snp-(nr_snp/2)+ii]-snps_pos[i_snp-(nr_snp/2) + ii-1]));
		}
		if( forwardArray[ii]+backwardArray[ii] +rprob > best_prob){
		  best_i = ii;
		  best_prob = forwardArray[ii]+backwardArray[ii]+rprob;
		}
	      }
	      /*new breakpoint location*/
	      fprintf(stderr,"%d %d -> %d ([%d %d] -> [%d %d]), prob (%f)\n",i_hap, i_snp , i_snp-(nr_snp/2)+best_i, left[0], left[1], right[0], right[1]  ,  best_prob);
	      if(best_i < (nr_snp/2)){
		for(int ii = i_snp-(nr_snp/2)+best_i;ii<i_snp;ii++){
		  admx_samples_inf[i_hap][ii] = (char)(10*right[0] +right[1]);
                }
	      }else{
		for(int ii = i_snp; ii< i_snp-(nr_snp/2)+best_i; ii++){
		  admx_samples_inf[i_hap][ii] = (char)(10*left[0] +left[1]);
		}
	      }
	    }
	  }
	}




#if 0
	for(int i = 0 ; i< nr_admx_samples; i++){
	  for(int j =1; j < nr_total_snps;j++){
	    if(admx_samples_inf[i][j] != admx_samples_inf[i][j-1]){
	      int jj=j+1;
	      while(admx_samples_inf[i][jj] == admx_samples_inf[i][jj-1] &&jj<nr_total_snps) 
		jj++;
	      if(jj < nr_total_snps && (jj-j) <= (nr_snp-window_padding)){
		for(int jjj =j; jjj<jj;jjj++){
		  admx_samples_inf[i][jjj] = admx_samples_inf[i][j-1];    
		}
	      }
	    }
	  }
	}
#endif
	for(int i = 0 ; i< nr_admx_samples; i++){
	  for(int j =1; j < nr_total_snps;j++){
	    if(admx_samples_inf[i][j] != admx_samples_inf[i][j-1]){
	      fprintf(fout," %02d:%d",admx_samples_inf[i][j-1],j);
	      // fprintf(fout,"Ind %d SNP %d %d trans %d to %d\n",i,j-1,snps_pos[j-1],admx_samples_inf[i][j-1],admx_samples_inf[i][j]);
	    }
	  }
	  fprintf(fout," %02d:%d",admx_samples_inf[i][nr_total_snps-1],nr_total_snps);
	  fprintf(fout,"\n");
	}
	fclose(fout);
	
	return 0;
	
}
char* getgeno( char* hap1, char* hap2, int nr_snp){
	char *gen =  (char*)calloc (nr_total_snps+1,sizeof(char));
	for(int i =0; i < nr_snp; i++){
		gen[i] = (char)( (hap1[i] -48)+ (hap2[i]-48) +48);
	}
}


void fix_ancestry_hap(int right_snp){
  for(int i =0; i < nr_admx_samples;i++){
    for(int j =starting_snp; j < right_snp;j++){
      double max = 0;
      int a1=0;
      int a2=0;
      for(int ii =0; ii < nr_anc;ii++){
	for(int jj =0; jj < nr_anc; jj++){
	  if(admx_inf_array[i][j][ii][jj] > max){
	    a1=ii;
	    a2=jj;
	    max = admx_inf_array[i][j][ii][jj];
	  }
	}
      }
      if(a1<a2){
	admx_samples_inf[i][j] = (char)(10*a1+a2);
      }else{
	admx_samples_inf[i][j] = (char)(10*a2+a1);
      }
    }
  }
}


void get_tot_acc(int right_snp, FILE *f){
	int tot_c=0;
	int tot_e=0; 
	//get accuracy
	for(int i = 0 ; i< nr_admx_samples; i++){
		
		for(int j =1; j < right_snp;j++){
			if(admx_samples_inf[i][j] != admx_samples_inf[i][j-1]){
				int jj=j+1;
				while(admx_samples_inf[i][jj] == admx_samples_inf[i][jj-1] &&jj<right_snp) 
					jj++;
				if(jj < right_snp && (jj-j) <= (nr_snp-window_padding)){
					for(int jjj =j; jjj<jj;jjj++){
						admx_samples_inf[i][jjj] = admx_samples_inf[i][j-1];    
					}
				}
			}
		}

		int c =0;
		int e=0;
		for(int j = starting_snp; j < right_snp;j++){
			char true_anc;
			int h1,h2;
			if(i == 0){
				h1 = 0;
				h2 = 1;
			}
			if( i == 1){
				h1 = 2;
				h2 = 3;
			}
			if( i == 2){
				h1 = 0;
				h2 = 2;
				
			}
				
			true_anc = (char)((admx_samples_anc[h1][j]-48)*10+(admx_samples_anc[h2][j]-48));
			if(admx_samples_anc[h1][j] < admx_samples_anc[h2][j]){
				true_anc = (char)((admx_samples_anc[h1][j]-48)*10+(admx_samples_anc[h2][j]-48));
			}else{
				true_anc = (char)((admx_samples_anc[h2][j]-48)*10+(admx_samples_anc[h1][j]-48));
			}
			
			//fprintf(stderr,"%d %d\n",true_anc,admx_samples_inf[i][j]); 
			if( true_anc != admx_samples_inf[i][j] ){
				e++;
			}else{
				c++;
			}
		}
		tot_c+=c;
		tot_e+=e;
		
		fprintf(stderr,"TOT: %d %d %d %f (%f)\n",i,e,e+c,(double)c/(e+c),(double)tot_c/(tot_e+tot_c));
#if 0				
		for(int j =1; j < right_snp;j++){
			bool true_bkpt = false;
			
			
			if(i<nr_admx_samples){
				if(admx_samples_anc[i][j] != admx_samples_anc[i][j-1]){
					true_bkpt =  true;
				}
			}else {
				if(admx_samples_anc[0][j] != admx_samples_anc[0][j-1]){
					true_bkpt =  true;
				}
				if(admx_samples_anc[1][j] != admx_samples_anc[1][j-1]){
					true_bkpt =  true;
				}
			}

				
				
				if(true_bkpt){
					/*find closest inferred bkpt*/
					
					int kleft = j;
					int cont = 1;
					while(cont && kleft > 0){
						if(admx_samples_inf[i][kleft] != admx_samples_inf[i][kleft-1]){
							cont = 0;
						}else{
							kleft--;
						}
					}
					int kright = j;
					cont = 1;
					while(cont && kright < nr_total_snps){
						if(admx_samples_inf[i][kright] != admx_samples_inf[i][kright-1] ){
							cont = 0;
						}else{
							kright++;
						}	
					}
					int pred_bkpt;
					int bp_dist;
					if((j-kleft) > (kright-j)){
						pred_bkpt = kright;
						if(pred_bkpt == nr_total_snps){
							pred_bkpt = nr_total_snps-1;
						}
						bp_dist = snps_pos[pred_bkpt] - snps_pos[j];    
					}else{
						pred_bkpt = kleft;
						bp_dist = snps_pos[j] -snps_pos[kleft];
					}
					
					
				//	fprintf(stdout,"%d %d %d %d %d\n",i,j,pred_bkpt,(int)fabs(j-pred_bkpt),bp_dist);
					
					
				}
		
		}
#endif		
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
  snps_map = (double*)calloc(nr_total_snps, sizeof(double));
  rewind(fin);	
  i=0;
  int tempint;
  char tempc;
  double tempdouble;
  //while (fscanf(fin,"%s\t%d\t%f\t%d\t%c\t%c\n",&myhap,&tempint,&(snps_map[i]),&(snps_pos[i]),&tempc,&tempc) > 0 ) {
  while (fscanf(fin,"%d\n",&(snps_pos[i])) > 0 ) {
    // fprintf(stderr,"%d\n",i);
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
  char myhap [MAX_SNPS];

  int i = 0;

  if(!fin){
    fprintf(stderr, "Could not open file: %s !!\n",filename);
    exit(1);
  }

  nr_admx_samples = 0;
  while (fgets(myhap,MAX_SNPS,fin)) {
    nr_admx_samples++;
  }

  admx_samples = (char**)calloc(nr_admx_samples+1, sizeof(char*));
  admx_samples_inf = (char**)calloc(nr_admx_samples+1 , sizeof(char*));
  admx_inf_array = (double****)calloc(nr_admx_samples+1,sizeof(double***));
  for(int i = 0; i < nr_admx_samples+1; i++)
    {
      admx_samples[i] = (char*)calloc((nr_total_snps+1) , sizeof(char));
      admx_samples_inf[i] = (char*)calloc((nr_total_snps+1) , sizeof(char));
      admx_inf_array[i] = (double***)calloc(nr_total_snps,sizeof(double**));
      for(int ii = 0; ii < nr_total_snps;ii++){
	admx_inf_array[i][ii] = (double**)calloc(nr_anc,sizeof(double*));
	for(int jj = 0; jj< nr_anc;jj++){
	  admx_inf_array[i][ii][jj] = (double*)calloc(nr_anc,sizeof(double));
	}
      }
    }
  rewind(fin);
  i=0;
  //while (fgets(admx_samples[i],MAX_SNPS,fin)) {                                                                                                     
  while (fgets(myhap,MAX_SNPS,fin)) {
    strncpy(admx_samples[i] , myhap, nr_total_snps*sizeof(char));
    i++;
  }
  fclose(fin);
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
    //fprintf(stderr,"%d\n", anc_map[ancestry].nr_hap);                                                                                                     
  }
  fprintf(stderr,"reading %d haps for ancestry %d\n",anc_map[ancestry].nr_hap,ancestry);
  anc_map[ancestry].hap_array = (char**)calloc(anc_map[ancestry].nr_hap , sizeof(char*));
  for(int i = 0; i < anc_map[ancestry].nr_hap; i++)
    {
      //fprintf(stderr,"%d\n",i);                                                                                                                           
      anc_map[ancestry].hap_array[i] = (char*)malloc((nr_total_snps+1) * sizeof(char));
    }
  rewind(fin);
  i=0;
  //while (fgets(anc_map[ancestry].hap_array[i],MAX_SNPS,fin)) {                                                                                            
  while (fgets(myhap,MAX_SNPS,fin)) {
    strncpy(anc_map[ancestry].hap_array[i] , myhap, nr_total_snps*sizeof(char));
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
	forward_hap_array = (double**)malloc(nr_snp*sizeof(double*));
	for(int i =0; i < nr_snp; i++)
	{
		forward_hap_array[i] = (double*)malloc(nr_founders*sizeof(double));
	}	
	backward_hap_array = (double**)malloc(nr_snp*sizeof(double*));
	for(int i =0; i < nr_snp; i++)
	{
		backward_hap_array[i] = (double*)malloc(nr_founders*sizeof(double));
	}
	
	//Uno
	array_uno_precomp1 = (double**)malloc(nr_founders*sizeof(double*));
	for(int i =0; i< nr_founders; i++)
	{
		array_uno_precomp1[i] = (double*)malloc(nr_founders*sizeof(double));
	}
	
	totalProbUnoMatrix_forward = (double***)malloc(nr_snp*sizeof(double**)); 
	totalProbUnoMatrix_backward = (double***)malloc(nr_snp*sizeof(double**));
	totalProbUnoMatrix_unk = (double**)malloc(nr_snp*sizeof(double*));
	
	for(int i =0; i < nr_snp; i++)
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
	
}
void alloc_mem_hmm(HAP_HMM *hmm)
{	
	hmm->states = (STATE**)malloc(nr_snp*sizeof(STATE*));
	for(int i =0; i < nr_snp; i++)
	{
		hmm->states[i] = (STATE*)malloc(nr_founders*sizeof(STATE));
	}
	hmm->pi = (double*)malloc(nr_founders*sizeof(double));
	hmm->pi_end = (double*)malloc(nr_founders*sizeof(double));
	hmm->trans = (double***)malloc((nr_snp-1)*sizeof(double**));
	for(int i =0; i < nr_snp-1; i++)
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
	for(int i =0; i < nr_snp; i++)
		free(forward_hap_array[i]);
	free(forward_hap_array);
	
	for(int i =0; i < nr_snp; i++)
		free(backward_hap_array[i]);
	free(backward_hap_array);
	//Uno
	for(int i =0; i< nr_founders; i++)
		free(array_uno_precomp1[i]);
	free(array_uno_precomp1);
	
	for(int i =0; i < nr_snp; i++)
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
	for(int i =0; i < nr_snp; i++)
		free(hmm->states[i]);
	free(hmm->states);
	free(hmm->pi);
	free(hmm->pi_end);
	for(int i =0; i < nr_snp-1; i++)
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
			fprintf(stderr,"ERROR forward_hap hap[%d]=%c not found\n\n", start_pos+j,hap[start_pos+j] );
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
				fprintf(stderr,"ERROR forward_hap hap[%d]=%c not found\n\n", start_pos+j,hap[start_pos+j] );
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
			      double probEM = 1;
			      if(anc_map[h].hap_array[index_hap][start_pos+ t+1] == '0' || anc_map[h].hap_array[index_hap][start_pos+ t+1] == '1'){
				probEM = anc_map[h].hmm.states[t+1][j].prob_em[ anc_map[h].hap_array[index_hap][start_pos+ t+1]-48 ];
			      }else{
				if(anc_map[h].hap_array[index_hap][start_pos+ t+1] == '0'){
				  probEM = 1;
				}
			      }
			      if(forward_hap_array[t][i] *
                                 anc_map[h].hmm.trans[t][i][j] *
                                 probEM*
				 backward_hap_array[t+1][j] / prob_of_hap < MIN_NR)
				{
				  fprintf(stderr,"ERROWRRRRRRRR!!! %e %e\n",MIN_NR, forward_hap_array[t][i] *
                                          anc_map[h].hmm.trans[t][i][j] *
                                          probEM*
					  backward_hap_array[t+1][j] / prob_of_hap);
				  exit(0);
				}
                              trans_temp[t][i][j] += forward_hap_array[t][i] *
                                anc_map[h].hmm.trans[t][i][j] *
                                probEM* //anc_map[h].hmm.states[t+1][j].prob_em[ anc_map[h].hap_array[index_hap][start_pos+t+1]-48 ]*                                          
                                backward_hap_array[t+1][j] / prob_of_hap ;
			    }
			}
		    }
		  for(int t = 0 ; t < nr_snp; t++){
		    if(anc_map[h].hap_array[index_hap][start_pos+t] == '0' || anc_map[h].hap_array[index_hap][start_pos+t]=='1'){
		      for(int i =0; i < nr_founders; i++)
			states_temp[t][i].prob_em[ anc_map[h].hap_array[index_hap][start_pos+t]-48 ] += forward_hap_array[t][i] * backward_hap_array[t][i] / prob_of_hap;
		    }
		  }
		  if(anc_map[h].hap_array[index_hap][start_pos+0] =='0' || anc_map[h].hap_array[index_hap][start_pos+0] == '1'){
		    for(int i =0; i < nr_founders; i++)
		      pi_temp[i] += anc_map[h].hmm.states[0][i].prob_em[ anc_map[h].hap_array[index_hap][start_pos+0]-48 ]*backward_hap_array[0][i] / prob_of_hap;
		  }
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
