#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

static int* append_edge(int* edges, int* nedge, int i, int j) {
    int buffer = 1024;
    int n = *nedge;
    if (n == 0) {
        edges = (int*)malloc(sizeof(int) * buffer * 2);
    }

    if ((n + 1) % buffer == 0) {
        edges = (int*)realloc(edges, sizeof(int)*((n+1)/buffer+1)*buffer*2);
    }

    edges[2*n+0] = i;
    edges[2*n+1] = j;
    n++;
    *nedge = n;

    return edges;
}

int* read_edgelist(char* filename, int* nnode, int* nedge) {
    FILE* fp = fopen(filename,"r");
    if(fp==NULL) {
        printf("Error!\n");
        exit(1);
    }

    int* edges = NULL;

    int i,j;
    char data[128];
    int nnode_temp = 0;
    *nedge = 0;
    while(fscanf(fp,"%d %d %s",&i,&j,data)==3) {
        if(nnode_temp<i) {
            nnode_temp=i;
        } else if(nnode_temp<j) {
            nnode_temp=j;
        }
        
        edges = append_edge(edges,nedge,i,j);
    }

    *nnode = nnode_temp+1;
    fclose(fp);

    return edges;
}

int* kinship_graphs_generator(int* nnode, int* nedge, int N, double rho, gsl_rng* rng) {
    int* edges  = NULL;
    int nc_max  = 100;
    int* groups = (int*)malloc(sizeof(int)*N*nc_max);
    int* g_list = (int*)malloc(sizeof(int)*2*N);
    int* n_list = (int*)malloc(sizeof(int)*N);

    for(int i=0; i<N; i++) n_list[i]=0;

    int redo_generator=1;
    while(redo_generator) {
    *nedge=0;
    if(edges!=NULL) {
        free(edges);
        edges=NULL;
    }
    for(int i_group=0; i_group<N; i_group++) {
        for(int i=0; i<2; i++) {
            if(gsl_rng_uniform_pos(rng)<rho) {
                int j_group = (int)(gsl_rng_uniform_pos(rng)*N);
                groups[j_group*nc_max+n_list[j_group]]=i_group*2+i;
                g_list[i_group*2+i]=j_group;
                n_list[j_group]++;
            } else {
                groups[i_group*nc_max+n_list[i_group]]=i_group*2+i;
                g_list[i_group*2+i]=i_group;
                n_list[i_group]++;
            }
        }
    }

    for(int i_group=0; i_group<N; i_group++) {
        for(int i=0; i<n_list[i_group]; i++) {
            for(int j=i+1; j<n_list[i_group]; j++) {
                int node_i = groups[i_group*nc_max+i];
                int node_j = groups[i_group*nc_max+j];

                edges = append_edge(edges,nedge,node_i,node_j);
            }
        }
    }

    redo_generator=0;
    for(int i_group=0; i_group<N; i_group++) {
        int node_i = i_group*2;
        int node_j = (int)(gsl_rng_uniform_pos(rng)*N)*2+1;
        int assigned=0;
        for(int i=0; i<N && assigned==0; i++) {
            if(g_list[node_i]==g_list[node_j] || g_list[node_j]==-1) {
                assigned=0;
            } else {
                assigned=1;
            }
            node_j = (node_j+2)%(2*N);
        }

        if(assigned){
            edges = append_edge(edges,nedge,node_i,node_j);
            g_list[node_i]=-1;
            g_list[node_j]=-1;
        } else {
            redo_generator=1;
            break;
        }
    }
    }

    *nnode=2*N;

    free(groups);
    free(g_list);
    free(n_list);

    return edges;
}

#if 0
int main() {
    int N=10000;
    double rho=1.0;
    int nsample=10000;

    int nnode,nedge;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,33798743);

    double k1_mean=0;
    double k2_mean=0;

    for(int i=0; i<nsample;i++) {
        int* edges  =  kinship_graphs_generator(&nnode, &nedge, N, rho, rng);
        int* degree = (int*)malloc(sizeof(int)*nnode);
        for(int j=0;j<nnode;j++) degree[j]=0;

        for(int j=0;j<nedge;j++) {
            degree[edges[j*2+0]]++;
            degree[edges[j*2+1]]++;
        }

        double k1=0;
        double k2=0;
        for(int j=0;j<nnode;j++) {
            k1+=degree[j];
            k2+=degree[j]*degree[j];
        }
        k1_mean += k1/nnode;
        k2_mean += k2/nnode;

        free(degree);
        free(edges);
    }
    k1_mean = k1_mean/nsample;
    k2_mean = k2_mean/nsample;

    printf("k1=%.8f k2=%.8f k1/k2=%.8f\n",k1_mean,k2_mean,k1_mean/k2_mean);

    gsl_rng_free(rng);
    return 0;
}
#endif
