#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "graph.h"
#include "time_evolution.h"

static void initial_state(int* sigma, int nnode, double p, int n, int type, gsl_rng* rng) {
    if(type==0) {
        for(int i=0;i<nnode;i++) {
            sigma[i]=0;
            if(gsl_rng_uniform_pos(rng)<p) sigma[i]=1;
        }
    } else if(type==1) {
        for(int i=0;i<nnode;i++) {
            sigma[i]=0;
            if(i<n) sigma[i]=1;
        }
    } else if(type==2) {
        for(int i=0;i<nnode;i++) {
            sigma[i]=0;
        }
        int dis=(int)(gsl_rng_uniform_pos(rng)*nnode);
        sigma[dis]=1;
    }
}

static int final_state(int nnode, int* sigma, double p) {
    int check=0;
    int n=0;
    for(int i=0;i<nnode;i++) {
        n+=sigma[i];
    }

    if(nnode*p<n) check=1;

    return check;
}

void print_state(int* sigma, int nnode) {
    for(int i=0;i<nnode;i++) {
        printf("%d ",sigma[i]);
    }
    printf("\n");
}

void save_state(FILE* file, int* sigma, int nnode) {
    for(int i=0;i<nnode;i++) {
        fprintf(file,"%d ",sigma[i]);
    }
    fprintf(file,"\n");
}

static double order_parameter(int* sigma, int nnode) {
    double I=0;
    for(int i=0;i<nnode;i++) {
        I+=sigma[i];
    }
    I = I/nnode;

    return I;
}

int main(int argc, char** argv) {
    double alpha=atof(argv[1]);
    double gamma=atof(argv[2]);
    double dt=atof(argv[3]);
    double T=atof(argv[4]);
    double p=atof(argv[5]);
    int nif = 10;
    int block_size=atoi(argv[6]);
    int nblock=atoi(argv[7]);
    unsigned long int seed=atoi(argv[8]);
    int nstep=(int)(T/dt);
    int conf_size=100;
    int nshow=nstep/conf_size;

    gsl_rng* rng=gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    int nnode,nedge;
    char filename[128] = "/hpc/home/jp549/src/ctQMC/C/projects/epidemic/network/test.edgelist";
    int* edges = read_edgelist(filename,&nnode,&nedge);


    double  infected_ratio=0;
    double* sigma_ave = (double*)malloc(sizeof(double)*nstep);
    double* temp_sigma = (double*)malloc(sizeof(double)*nstep);
    for(int i=0;i<nstep;i++) {
        sigma_ave[i]=0;
        temp_sigma[i]=0;
    }
    int* sigma = (int*)malloc(sizeof(int)*nnode);
    int* conf = (int*)malloc(sizeof(int)*nnode*(conf_size+1));
    for(int i=0;i<(nnode*(conf_size+1));i++) conf[i]=0;

    for(int k=0;k<nblock;k++){
        int n=0;
        time_t start_time = clock();
        FILE* file_conf = fopen("conf.txt","a");

        double ntrial_ave=0;
        int ntrial=0;
        for(int i_block=0;i_block<block_size;) {
            n=0;
            initial_state(sigma,nnode,p,nif,2,rng);

            int i_show=0;
            for(int i_node=0;i_node<nnode;i_node++)
                conf[i_show*nnode+i_node]=sigma[i_node];
            i_show++;

            infected_ratio = order_parameter(sigma,nnode);
            temp_sigma[n] = infected_ratio;
            n++;
            for(int i=0;i<nstep;i++) {
                kernel_linear(alpha,gamma,dt,nnode,nedge,sigma,edges,rng);
                if((i+1)%nshow==0) {
                    infected_ratio = order_parameter(sigma,nnode);
                    temp_sigma[n] = infected_ratio;
                    n++;

                    for(int i_node=0;i_node<nnode;i_node++)
                        conf[i_show*nnode+i_node]=sigma[i_node];

                    i_show++;
                }
            }
            
            ntrial++;
            if(final_state(nnode,sigma,p)) {
                for(int j=0;j<n;j++) sigma_ave[j]+=temp_sigma[j];

                ninfection_count_plus_one();
                nrecover_count_plus_one();

                printf("i_block = %d | trial = %d \n",i_block,ntrial);
                ntrial_ave+=ntrial;
                ntrial=0;

                for(i_show=0;i_show<(conf_size+1);i_show++)
                    save_state(file_conf,&(conf[i_show*nnode]),nnode);

                i_block++;
            }

        }

        time_t end_time = clock();
        FILE* file_t = fopen("times.txt","w");
        FILE* file_s = fopen("series.txt","a");
        FILE* file_g = fopen("global.txt","a");
        printf("-----------------------\n");
        printf(" t    |    I/N \n");
        for(int i=0;i<n;i++) {
            sigma_ave[i] = sigma_ave[i]/block_size;
            printf("%.4lf  %.12lf\n",nshow*dt*i,sigma_ave[i]);
            fprintf(file_t,"%.4lf ",nshow*dt*i);
            fprintf(file_s,"%.12e ",sigma_ave[i]);
            sigma_ave[i] = 0;
        }
        fprintf(file_t,"\n");
        fprintf(file_s,"\n");
        printf("time for this block = %.2lf(sec)\n",(double)(end_time-start_time)/CLOCKS_PER_SEC);
        double ninfection = ninfection_ave_value();
        double nrecover = nrecover_ave_value();
        ntrial_ave = ntrial_ave/block_size;
        fprintf(file_g,"%.12e %.12e %.12e\n",ninfection,nrecover,ntrial_ave);

        print_ninfection();
        print_nrecover();

        //fclose(file_conf);
        fclose(file_t);
        fclose(file_s);
        fclose(file_g);
    }

    free(edges);
    free(sigma);
}
