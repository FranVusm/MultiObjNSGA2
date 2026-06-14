/* Routine for evaluating population members  */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "global.h"
# include "rand.h"


int isADepot(int name, problem_instance pi){
    int i, j;
    int count=0;
    for (i=0; i<pi.nDepots; i++){
        for (j=0; j<3; j++){
            if(name == pi.depots[i].names[j])
                return 1; /*Es un depot*/
        }
        count++;
    }
    for (i=0; i<pi.nClients; i++){
        if(name == pi.clients[i].name)
            return 0; /*No es un depot*/
        count++;
    }
    return -1;
}

int getNodeId(int name, problem_instance pi){
    int i, j;
    int count=0;
    for (i=0; i<pi.nDepots; i++){
        for (j=0; j<3; j++){
            if(name == pi.depots[i].names[j])
                return count;
        }
        count++;
    }
    for (i=0; i<pi.nClients; i++){
        if(name == pi.clients[i].name)
            return count;
        count++;
    }
    return -1;
}

int getNodeIndex(int name, problem_instance pi){
    int i, j;
    int count=0;
    for (i=0; i<pi.nDepots; i++){
        for (j=0; j<3; j++){
            if(name == pi.depots[i].names[j])
                return count;
        }
        count++;
    }
    count=0; /*Se cuenta desde cero otra vez*/
    for (i=0; i<pi.nClients; i++){
        if(name == pi.clients[i].name)
            return count;
        count++;
    }
    return -1;
}



void printProblemInstance(problem_instance *pi){
    int i, j;
    printf("Depots: %d\n", pi->nDepots);
    for (i=0; i<pi->nDepots; i++){
        printf("%d: ", pi->depots[i].id);
        for (j=0; j<3; j++){
            printf("%ld ", pi->depots[i].names[j]);
        }
        printf("\n");
    }

    printf("Clients:\n");
    for (i=0; i<pi->nClients; i++){
        printf("%d , %d, %d\n", pi->clients[i].id, pi->clients[i].name, pi->clients[i].demand);
    }

    printf("Trucks: %d\n", pi->nTrucks);
    for (i=0; i<pi->nTrucks; i++){
        printf("Truck: %d, Capacity: %d \n", pi->trucks[i].id, pi->trucks[i].capacity);
        printf("dStart: ");
        for (j=0; j<pi->trucks[i].ndStart; j++){
            printf("%d ", pi->trucks[i].dStart[j]);
        }
        printf("dIn: ");
        for (j=0; j<pi->trucks[i].ndIn; j++){
            printf("%d ", pi->trucks[i].dIn[j]);
        }
        printf("dEnd: ");
        for (j=0; j<pi->trucks[i].ndEnd; j++){
            printf("%d ", pi->trucks[i].dEnd[j]);
        }

        printf("\n");
    }


    return;

}
