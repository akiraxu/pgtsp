#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <math.h>
#include <random>
#include <time.h>
#include <unordered_set>
#include "mpi.h"

using namespace std;

int paramN; //# of sites
int *paramD; // Distance matrix
int paramM;
int paramK; //Max round
int paramR; //No change for
int paramS; //# of survivors
int paramU; //# of mutated individuals
int paramZ; // mutated pairs
default_random_engine *gen;

void initEverything(int n, int m, int k, int r, float s, float u, float z){
    srand(time(NULL));
    gen = new default_random_engine(rand());
    paramN = n;
    paramM = m;
    paramK = k;
    paramR = r;
    paramS = s * paramM;
    paramU = u * paramM;
    paramZ = z * paramN;
}

void setParams(int *arr){
    srand(time(NULL));
    gen = new default_random_engine(rand());
    paramN = arr[0];
    paramM = arr[1];
    paramK = arr[2];
    paramR = arr[3];
    paramS = arr[4];
    paramU = arr[5];
    paramZ = arr[6];
}

int randInt(){
    uniform_int_distribution<> dis(0, 2147483647);
    return dis(*gen);
}

int twoCityDistance(int a, int b){
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
    for(int i = paramS * paramN; i < paramN * paramM; i += paramN){
        if(flag){
            x = randInt() % paramS;
            y = randInt() % paramS;
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
        mutation(paths + (randInt() % paramM) * paramN);
    }
}

void selection(int *distances, int *survivors){
    unsigned long long int totalDistance = 0;
    unsigned long long int tries = 0;
    for(int i = 0; i < paramM; i++){
        totalDistance += distances[i];
    }
    uniform_int_distribution<unsigned long long int> dis(0, totalDistance - 1);
    unordered_set<int> used;
    for(int i = 0; i < paramS; i++){
        int n = -1;
        while(true){
            tries++;
            n = -1;
            unsigned long long int k = dis(*gen);
            unsigned long long int t = 0;
            for(int j = 0; j < paramM; j++){
                if(t <= k && t + distances[j] > k){
                    n = j;
                    break;
                }
                t += distances[j];
            }

            if(used.count(n) == 0){
                break;
            } 
        }
        survivors[i] = n;
        used.insert(n);
    }
    cout << "avg try per survivor: " << 1.0 * tries / paramS << endl;
}

void copySelection(int *from, int *to, int *survivors){
    for(int i = 0; i < paramS; i++){
        for(int j = 0; j < paramN; j++){
            to[i * paramN + j] = from[survivors[i] * paramN + j];
        }
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
void outPathToFile(int *arr, string fn){
    ofstream out(fn);
    for(int i = 0; i < paramM; i++){
        for(int j = 0; j < paramN; j++){
            out << arr[i * paramN + j] << " ";
        }
        out << endl;
    }
    out.close();
}

int minDistanceIndex(int *arr){
    int index = 0;
    int min = arr[0];
    for(int i = 0; i < paramM; i++){
        if(arr[i] < min){
            min = arr[i];
            index = i;
        }
    }
    return index;
}

int main(int argc,char* argv[]){

    int id;
	int p;
	double wtime;

    initEverything(stoi(argv[1]), stoi(argv[2]), stoi(argv[3]), stoi(argv[4]), stof(argv[5]), stof(argv[6]), stof(argv[7]));
    paramD = new int[paramN * paramN];

	MPI::Init(argc, argv); //  Initialize MPI.
	p = MPI::COMM_WORLD.Get_size(); //  Get the number of processes.
	id = MPI::COMM_WORLD.Get_rank(); //  Get the individual process ID.
    
    if(id = 0){
        string fn(argv[8]);
        ifstream in(fn);
        int n = 0;
        while(!in.eof())
        {
            in >> paramD[n];
            n++;
        }
        if(n != paramN * paramN + 1){
            cout << "Wrong distance matrix! " << n << " Exit now..." << endl;
            exit(-1);
        }
    }

    MPI_Bcast(paramD,paramN * paramN,MPI_INT,0,MPI_COMM_WORLD);

    int *paths = new int[paramN * paramM];
    int *paths2 = new int[paramN * paramM];
    int *distances = new int[paramM];
    int *survivors = new int[paramM];

    initPopulation(paths);

    for(int i = 0; i < paramN * paramM; i += paramN){
        randPath(paths + i);
    }

    int min = -1;

    int * minPath = new int[paramN];

    int countR = 0;

    for(int whatever = 0; whatever < paramK; whatever++){
        for(int i = 0; i < paramM; i++){
            distances[i] = roundDistance(paths + i * paramN);
        }
        
        int currmin = distances[minDistanceIndex(paths)];

        if(min == -1 || currmin < min){
            min = currmin;
            countR = 0;
            for(int i = 0; i < paramN; i++){
                minPath[i] = paths[currmin * paramN + i];
            }
        }else{
            countR++; 
        }

        if(countR > paramR){
            break;
        }

        cout << "Current round min distance: " << currmin << ", all time min: " << min << ", keep for " << countR << " rounds." << endl;

        selection(distances, survivors);

        int *tempPtr = paths;
        paths = paths2;
        paths2 = tempPtr;

        copySelection(paths2, paths, survivors);

        fillCrossover(paths);

        doMutation(paths);
    }
    
    cout << "min path: ";
    for(int i = 0; i < paramN; i++){
        cout << minPath[i] << " ";
    }
    cout << endl;

    MPI::Finalize();
}
