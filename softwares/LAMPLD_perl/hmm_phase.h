#include <stdio.h>
#define EPS 2.71828


/*DEFAULTS most not used*/
const double MIN_NR   = 1e-100;// used to compare probabilities
const double TRANSITION_PSEUDOCOUNT   = 0.00006;	
const double EMISSION_PSEUDOCOUNT     = 0.001;	
const double MIN_LOG  =  -10000000;
const double MIN_PROB = 0.0000000000001; 										//never let any prob get lower
const double MAJOR_PROB = 0.9;
const int MAX_LINE = 1000;			// length of sequences
const int TRUE = 1;               	// flag
const int FALSE = 0;              	// flag
/*GLOBALS USED in defining the HMM structure*/
int nr_snp; 	 					/*number of snp's*/
int nr_gen; 	 					/*number of genotypes*/
int nr_founders; 					/*number of founders*/
char** gen_array;					//array of genotypes
/*STRUCTS*/
typedef struct state
{
  double prob_em[3]; 				//prob of emitting 0 or 1 (or 2 in case of labeling hmm)
  int major; 						//takes values 0 or 1 indicating which alele is major
} STATE; 

/*
  This struct defines the HMM constructed as a sequence of nr_snp layers, each layer containing nr_founder states;
  there exists a transition from state i in layer l to state i' in layer l' iff l'= l+1	
*/
typedef struct hap_hmm
{
  STATE**	states;	        		//array of states nr_founders * nr_snp 
  STATE 	start;					//start state
  STATE 	end;						//end state
  double* pi;						//prob of transition from start to state i in first layer (size = nr_founders)
  double*** trans;				//prob of transition size = nr_snps*nr_founders*nr_founders;
  //trans[i][j][k] = prob of transition from state j in layer i to state k in layer i+1
  //trans[i][j][k] = t(states[i][j],states[i+1][k]);
  double* pi_end;					//prob of transition from state i in last layer to end_state
}HAP_HMM;	

typedef struct ANC_HMM //Used to describe the id of an ancestor, and it's corresponding HMM
{
  char name[20];
  int nr_hap;		//Count of haplotypes used to train HMM
  char** hap_array; //Block-based hap_array
  HAP_HMM hmm;
  double **forward_hap_array;
  double **backward_hap_array;
}ANC_HMM;



/*FUNCTIONS*/
/* function that allocates memory to given hmm
	WARNING!!! set nr_snp,nr_gen,nr_founders first
*/

void alloc_mem_hmm(HAP_HMM *hmm); 
/* 
function that initializez hmm with default values
*/

void def_val_hmm(HAP_HMM *hmm); 
/* 
	function that dealocates mem for hmm 
*/
void free_mem_hmm(HAP_HMM *hmm); 
void print_hmm(HAP_HMM hmm);
void print_hmm_file(HAP_HMM hmm, FILE *fout);
int read_genotypes (char*** gen_array, int* nr_snp, int* nr_gen);
/* prob of generating gen_0,gen_1 ...gen_i-1 and ending at state a,b at level i*/
double forward(int a, int b, int i , int gen ,HAP_HMM hmm);
/* prob of starting at a,b  and generating gen_i+1,gen_i+2 ...gen_nr_snp-1 */
double backward(int a, int b, int i , int gen ,HAP_HMM hmm);
