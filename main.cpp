#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "utils.h"
#include "HMMProblem.h"
using namespace std;

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void createByRef(double** &val, int sz1, int sz2);
void destroyByRef(double** &val, int sz1);

int main (int argc, char ** argv) {
	clock_t tm0 = clock();

//	//
//	// Test repeated memory grabbing and freeing
//	//
//	int n1 = 10000, n2 = 20000, rep = 1, i,j,k;
//	double **lotsmem = NULL;
//	for(i=0; i<rep; i++) {
////		lotsmem = Malloc(double*, n1);
//        lotsmem = new double* [n1];
//		for(j=0; j<n1; j++) {
////			lotsmem[j] = Malloc(double, n2);
//			lotsmem[j] = new double [n2];
//			for(k=0; k<n2; k++)
//				lotsmem[j][k] = 2.123;
//		}
//		for(k=0; k<1000000; k++) ;
//		for(j=0; j<n1; j++) {
////			free(lotsmem[j]);
//			delete [] lotsmem[j];
//            lotsmem[j] = NULL;
//        }
////		free(lotsmem);
//		delete [] lotsmem;
//        lotsmem = NULL;
//	}
//	// Max mem grabbed for n1=10000, n2=10000, rep=1, is 40MB

//	//
//    // substr check for awk
//    //
//    string s = "We think in generalities, but we live in details.";
//    if(s.substr(s.length()-1,1)==".")
//        printf("ooh\n");

//    //
//    // array if random numbers and percentiles test
//    //
//    int folds = 3;
//    int n = 100;
//    int ar[n];
//    int cnt[folds];
//    for(int i=0; i<folds; i++) cnt[i] = 0;
//    srand ( time(NULL) );
//    for(int i=0; i<n; i++) {
//        ar[i] = rand() % folds;
//        cnt[ ar[i] ]++;
//    }
//    for(int i=0; i<folds; i++)
//        printf("%d-th fold contains %d elements\n",i+1,cnt[i]);
    
//    //
//    // parsing multiple kcs
//    //
//    char str[] = "Using simple numbers-1~~~Using large numbers-1~~Write expression, negative slope-1";
//    char* kc;
//    kc  = strtok(str,"~\n\r");
//    while(kc != NULL) {
//        fprintf(stdout,"%s\n",kc);
//        kc  = strtok(NULL,"~\n\r");
//    }

//    //
//    // tracing alpha and beta
//    //
//    struct data** dt = Calloc(struct data*, 1);
//    dt[0] = Calloc(struct data, 1);
//    NPAR bigo01[42] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo02[29] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0};
//    NPAR bigo03[32] =  {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo04[1] =  {0};
//    NPAR bigo05[5] =  {0,0,0,0,0};
//    NPAR bigo06[8] =  {1,0,0,0,1,0,0,0};
//    NPAR bigo07[5] =  {0,0,0,0,0};
//    NPAR bigo08[38] =  {0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0};
//    NPAR bigo09[33] =  {0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo10[1] =  {0};
//    NPAR bigo11[26] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0};
//    NPAR bigo12[1] =  {1};
//    NPAR bigo13[15] =  {1,0,0,0,1,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo14[27] =  {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0};
//    NPAR bigo15[6] =  {1,0,1,0,0,0};
//    NPAR bigo16[30] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo17[34] =  {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0};
//    NPAR bigo18[1] =  {0};
//    NPAR bigo19[11] =  {1,0,0,0,0,0,0,1,0,0,1};
//    NPAR bigo20[34] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0};
//    NPAR bigo21[34] =  {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo22[15] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo23[28] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo24[7] =  {0,1,0,0,0,0,0};
//    NPAR bigo25[28] =  {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0};
//    NPAR bigo26[8] =  {0,0,0,0,0,1,0,0};
//    NPAR bigo27[23] =  {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1};
//    NPAR bigo28[41] =  {1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0};
//    NPAR bigo29[37] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo30[21] =  {0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0};
//    NPAR bigo31[42] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo32[29] =  {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
//    NPAR bigo33[36] =  {1,1,1,1,0,0,1,0,0,1,1,1,1,0,1,0,1,0,0,1,0,0,1,1,1,0,0,1,0,1,0,0,0,0,0,0};
//    NPAR bigo34[27] =  {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo35[36] =  {0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo36[46] =  {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0};
//    NPAR bigo37[3] =  {1,1,1};
//    NPAR bigo38[6] =  {0,0,0,0,0,0};
//    NPAR bigo39[49] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0};
//    NPAR bigo40[28] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
//    NPAR bigo41[27] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo42[21] =  {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo43[33] =  {0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    NPAR bigo44[34] =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};
//    NPAR bigo45[5] =  {0,0,1,0,0};
//    NPAR *bigo[45] = {bigo01,bigo02,bigo03,bigo04,bigo05,bigo06,bigo07,bigo08,bigo09,bigo10,bigo11,bigo12,bigo13,bigo14,bigo15,bigo16,bigo17,bigo18,bigo19,bigo20,bigo21,bigo22,bigo23,bigo24,bigo25,bigo26,bigo27,bigo28,bigo29,bigo30,bigo31,bigo32,bigo33,bigo34,bigo35,bigo36,bigo37,bigo38,bigo39,bigo40,bigo41,bigo42,bigo43,bigo44,bigo45};
//    int sizes[45] = {42,29,32,1,5,8,5,38,33,1,26,1,15,27,6,30,34,1,11,34,34,15,28,7,28,8,23,41,37,21,42,29,36,27,36,46,3,6,49,28,27,21,33,34,5};
//    NUMBER big_ll = 0.0, big_loglik = 0.0;
//    for(int ii=0; ii<45; ii++) {
//        dt[0]->obs =   bigo[ii];
//        dt[0]->n = sizes[ii];
//        dt[0]->cnt = 0;
//        dt[0]->alpha = NULL;
//        dt[0]->beta = NULL;
//        dt[0]->gamma = NULL;
//        dt[0]->xi = NULL;
//        dt[0]->p_O_param = 0.0;
//        struct param* param = Malloc(struct param, 1);
//        param->nS = 2;
//        param->nO = 2;
//        NUMBER *PI = init1DNumber(param->nS); PI[0] = 1.0; PI[1] = 0.0;
//        NUMBER **A = init2DNumber(param->nS,param->nS); A[0][0] = 1.0; A[0][1] = 0.0; A[1][0] = 1.0; A[1][1] = 0.0;
//        NUMBER **B = init2DNumber(param->nS,param->nO); B[0][0] = 0.916; B[0][1] = 0.084; B[1][0] = 0.3; B[1][1] = 0.7;
//
//        // compute pLo
//        HMMProblem::computeAlphaAndPOParam(1, dt, PI, A, B, param->nS);
//        printf("- Alpha -\n");
//    //    for(int t=0; t<dt[0]->n; t++)
//    //        for(int i=0; i<param->nS; i++)
//    //            printf("%8.6f%s",dt[0]->alpha[t][i],((i==0)?"\t":"\n"));
//        printf("LL:    %8.6f\n",-fsafelog(dt[0]->p_O_param));
//        big_loglik += -fsafelog(dt[0]->p_O_param);
//        // simulate
//        NUMBER *local_pred = init1DNumber(param->nO);
//        NUMBER *pL = init1DNumber(param->nS);
//        NUMBER *pLe = init1DNumber(param->nO);
//        NUMBER pLe_denom;
//        NUMBER ll = 0.0;
//        NUMBER prob;
//        for(int i=0;i<param->nO; i++) pL[i] = PI[i];
//        NPAR o,m,i,j, metrics_target_obs=0, isTarget;
//        for(int t=0; t<dt[0]->n; t++) {
//            o = dt[0]->obs[t];//[t];
//            isTarget = (metrics_target_obs == o);
//            for(m=0; m<param->nO; m++) local_pred[m] = 0.0;
//            // produce prediction and copy to result
//            for(m=0; m<param->nO; m++)
//                for(i=0; i<param->nS; i++)
//                    local_pred[m] += pL[i] * B[i][m];
//    //        for(m=0; m<param->nO; m++)
//    //            printf("%10.8f%s",local_pred[m],(m<(param->nO-1))?"\t":"\n");
//            // update p(L)
//            pLe_denom = 0.0;
//            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
//            for(i=0; i<param->nS; i++) pLe_denom += pL[i] * B[i][o];
//            for(i=0; i<param->nS; i++) pLe[i] = pL[i] * B[i][o] / safe0num(pLe_denom);
//            // 2. L = (pLe'*A)';
//            for(i=0; i<param->nS; i++) pL[i] = 0.0;
//            for(j=0; j<param->nS; j++)
//                for(i=0; i<param->nS; i++)
//                    pL[j] += pLe[i] * A[i][j];
//            prob = safenum(local_pred[metrics_target_obs]);
//            ll -= fsafelog(  prob)*   isTarget  +  fsafelog(1-prob)*(1-isTarget);
//        } // for all data
//        printf("LLvulg:%8.6f\n",ll);
//        printf("diff: %10.8f\n",(-fsafelog(dt[0]->p_O_param)-ll));
//        big_ll += ll;
//        RecycleFitData(1, dt, param);
//        free(PI);
//        free2DNumber(A, param->nS);
//        free2DNumber(B, param->nS);
//        free(param);
//        free(local_pred);
//        free(pLe);
//        free(pL);
//    }
//    free(dt[0]);
//    free(dt);
//    printf("big_ll %8.6f, big_ll %8.6f, diff=%8.6f\n",big_loglik, big_ll,big_loglik-big_ll);

    //
    // safe logit and sigmoid
    //
//    NUMBER res = 0;
//    for(int i=0; i<100000; i++) res = sigmoid(1.3); //0.0018
//    for(int i=0; i<100000; i++) res = 1 / ( 1 + exp(-1.3)); //0.0012
//    for(int i=0; i<100000; i++) res = sigmoid(2.3); //   (fast), 0.0012 (faster)
    
//    for(int i=0; i<100000; i++) res = logit(0.7); // 0.0037
//    for(int i=0; i<100000; i++) res = logit(0.7); // 0.0008 (fast) 0.0005 (faster)
//    for(int i=0; i<100000; i++) res = logit(0.7); // 0.0037
    
//    printf("logit(.4)=%8.6f, logit(1)=%8.6f, logit(.4)=%8.6f, logit(1)=%8.6f\n",logit(.4),logit(1),logit(.4),logit(1.0));

//    int n = 10000;
//    for(int i=0; i<=n; i++) {
////        printf("logit(%8.6f)=%8.6f, logit(%8.6f)=%8.6f\n",(NUMBER)i/n,logit((NUMBER)i/n),(NUMBER)i/n,logit((NUMBER)i/n));
//        printf("sigmoid(%8.6f)=%8.6f, sigmoid(%8.6f)=%8.6f\n",(NUMBER)i/n,sigmoid((NUMBER)i/n),(NUMBER)i/n,sigmoid((NUMBER)i/n));
//    }
//	for(int i=1;i<argc;i++)
//	{
//		if(argv[i][0] != '-') break; // end of options stop parsing
//		if(++i>=argc)
//			exit(1);
//		switch(argv[i-1][1])
//		{
//			case 'g':
//                char *ch;
//                int j=0;
//                int n = 0;
//                bool start = true;
//                ch = strtok(argv[i],",-\n\t\r");
//                while( ch != NULL) {
//                    if(start) n++;
//                    fprintf(stdout,"%s %d\n",(start)?"start":"finish",atoi(ch));
//                    ch = strtok(NULL,",-\n\t\r");
//                    j++;
//                    start = !start;
//                }
//                fprintf(stdout,"%d groups\n",n);
//                break;
//		}
//	}
    int sz1 = 3, sz2 = 4;
    double** val;
    createByRef(val, sz1, sz2);
    destroyByRef(val, sz1);
    
    
    printf("done in %8.6f seconds\n",(double)(clock()-tm0)/CLOCKS_PER_SEC);
    return 0;
}

void createByRef(double** &val, int sz1, int sz2) {
    val = Malloc(double*, sz1);
    for(int i=0; i<sz1; i++)
        val[i] = Malloc(double, sz2);
        
    for(int i=0; i<sz1; i++)
        for(int j=0; j<sz1; j++)
            val[i][j] = 1.0;
}

void destroyByRef(double** &val, int sz1) {
    for(int i=0; i<sz1; i++)
        free(val[i]);
    free(val);
    val = NULL;
}
