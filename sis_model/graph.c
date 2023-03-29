#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Appends an edge to the edges list, and resizes the memory if necessary.
 *
 * This function appends an edge (i, j) to the given edges list, reallocating memory in
 * chunks of buffer size if needed.
 *
 * @param edges Pointer to the integer array representing the edge list.
 * @param nedge Pointer to the integer representing the number of edges.
 * @param i The first node of the edge to be added.
 * @param j The second node of the edge to be added.
 * @return The updated edges list with the new edge appended.
 */
static int* append_edge(int* edges, int* nedge, int i, int j) {
    int buffer=1024;
    int n = *nedge;
    if(n==0) {
        edges = (int*)malloc(sizeof(int)*buffer*2);
    }

    if((n+1)%buffer==0) {
        int* edges_temp = (int*)malloc(sizeof(int)*((n+1)/buffer+1)*buffer*2);
        for(int k=0;k<2*n;k++) {
            edges_temp[k] = edges[k];
        }
        free(edges);
        edges = edges_temp;
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

#ifdef test_graph_c
int main() {
    int nnode,nedge;
    char filename[128] = "../../networks/jupyters/test.edgelist";
    int* edges = read_edgelist(filename,&nnode,&nedge);

    for(int i=0;i<nedge;i++) {
        printf("%d %d \n",edges[2*i],edges[2*i+1]);
    }

    printf("nnode=%d nedge=%d \n",nnode,nedge);

    free(edges);

    return 0;
}
#endif
