#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <math.h>

using namespace std;

int paramN = 100; //# of sites
int *paramD; // Distance matrix
int paramM = 100000;
int paramK = 1000; //Max round
int paramR = 10; //No change for
int paramS = 0.75 * paramM; //# of survivors
int paramU = 0.1 * paramM; //# of mutated individuals
int paramZ = 0.05 * paramN; // mutated pairs

int twoCityDistance(int a, int b,){
    return paramD[a * paramN + b];
}

int roundDistance(int *path){
    int total = 0;
    for(int i = 0; i < paramN; i++){
        total += twoCityDistance(path[i], path[i+1 == paramN ? 0 : i+1]);
    }
    return total;
}

void findCrossoverMissing(int *path, int *missing){
    int m = 0;
    for(int i = 0; i < paramN; i++){
        bool flag = false;
        for(int j = 0; j < paramN / 2; j++){
            if(path[j] == i){
                flag = true;
            }
        }
        if(!flag){
            missing[m] = i;
            m++;
        }
    }
}

void findCrossoverOrder(int *path, int *order){
    int o = 0;
    for(int i = 0; i < paramN; i++){
        for(int j = paramN / 2; j < paramN; j++){
            if(path[j] == i){
                order[o] = j;
                o++;
                break;
            }
        }
    }
}

void crossoverOne(int *x, int *y, int *v){
    int * missing = new int[paramN - paramN / 2];
    int * order = new int[paramN - paramN / 2];
    for(int i = 0; i < paramN / 2; i++){
        v[i] = x[i];
    }
    findCrossoverMissing(x, missing);
    findCrossoverOrder(y, order);
    for(int i = 0; i < paramN - paramN / 2; i++){
        v[order[i]] = missing[i];
    }
    delete[] missing;
    delete[] order;
}

void crossover(int *x, int *y, int *v, int *w){
    crossoverOne(x, y, v);
    crossoverOne(y, x, w);
}

void fillCrossover(int *paths){
    bool flag = true;
    int x, y;
    for(int i = (paramM - paramS) * paramN; i < paramN * paramM; i += paramN){
        if(flag){
            x = rand() % (paramM - paramS);
            y = rand() % (paramM - paramS);
            crossoverOne(paths + x * paramN, paths + y * paramN, paths + i);
        }else{
            crossoverOne(paths + y * paramN, paths + x * paramN, paths + i);
        }
        flag = !flag;
    }
}

void mutation(int *path){
    for(int i = 0; i < paramZ; i++){
        int a = rand() % paramN;
        int b = rand() % paramN;
        int n = path[a];
        path[a] = path[b];
        path[b] = n;
    }
}

void doMutation(int *paths){
    for(int i = 0; i < paramU; i++){
        mutation(paths + (rand() % paramM) * paramN);
    }
}

void selection(int *distances, int *survivors){
    int totalDistance = 0;
    for(int i = 0; i < paramM; i++){
        totalDistance += distances[i];
    }
    for(int i = 0; i < paramM - paramS; i++){
        int n = -1;
        while(true){
            n = -1;
            int k = rand() % totalDistance;
            int t = 0;
            for(int j = 0; j < paramN; j++){
                if(t < k && t + distances[j] >= k){
                    n = j;
                    break;
                }
                t += distances[j];
            }
            bool flag = true;
            for(int j = 0; j < i; j++){
                if(survivors[j] == n){
                    flag = false;
                }
            }
            if(flag){
                break;
            }
        }
        survivors[i] = n;
    }
}

void copySelection(int *from, int *to, int *survivors){
    for(int i = 0; i < paramM - paramS; i++){
        to[i * paramN] = from[survivors[i] * paramN];
    }
}

void randPath(int *path){
    int t, k;
    for(int i = 0; i < paramN; i++){
        t = path[i];
        k = rand() % paramN;
        path[i] = path[k];
        path[k] = t;
    }
}

void initPopulation(int *arr){
    for(int i = 0; i < paramN * paramM; i++){
        arr[i] = i % paramN;
    }
}

int main(int argc,char* argv[]){
    string fn(argv[1]);
    ifstream in(fn);
    paramD = new int[paramN * paramN];
    int n = 0;
    while(!File.eof())
    {
        File >> arr[n];
        n++;
    }
    if(n != paramN * paramN){
        cout << "Wrong distance matrix! Exit now..." << endl;
        exit(-1);
    }

    int *paths = new int[paramN * paramM];
    int *paths2 = new int[paramN * paramM];
    int *distances = new int[paramM];
    int *survivors = new int[paramM - paramS];

    initPopulation(paths);

    for(int i = 0; i < paramN * paramM; i += paramN){
        randPath(paths + i);
    }

    for(int i = 0; i < paramM; i++){
        distances[i] = roundDistance(paths + i * paramN);
    }

    selection(distances, survivors);

    int *tempPtr = paths;
    paths = paths2;
    path2 = tempPtr;

    copySelection(paths2, paths, survivors);

    fillCrossover(paths);

    doMutation(paths);

}