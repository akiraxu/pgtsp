#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <math.h>
#include <random>
#include <sys/time.h>
#include <unordered_set>
#include "mpi.h"

using namespace std;

int paramN; //# of sites
int *paramD; // Distance matrix
int paramM; //Population
int paramK; //Max round
int paramR; //No change for
int paramS; //# of survivors
int paramU; //# of mutated individuals
int paramZ; //Mutated pairs
int paramE; //Exchange once how many round
int paramP; //How many processes
int myID;

float ss, uu, zz;
default_random_engine *gen;

void initEverything(int n, int m, int k, int r, float s, float u, float z, int e){
    paramN = n;
    paramM = m;
    paramK = k;
    paramR = r;
    paramS = s * paramM;
    paramU = u * paramM;
    paramZ = z * paramN;
    paramE = e;

    ss = s;
    uu = u;
    zz = z;
}

void updateM(int m){
    paramM = m;
    paramS = ss * paramM;
    paramU = uu * paramM;
}

void initRand(int d){
    timeval t;
    gettimeofday(&t, NULL);
    srand(t.tv_usec % (time(NULL) + d));
    cout << d << ":" << rand() << endl;
    gen = new default_random_engine(rand());
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
    //cout << "avg try per survivor: " << 1.0 * tries / paramS << endl;
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

void randomExchangeList(int *list){
    unordered_set<int> used;
    for(int i = 0; i < paramM / (paramP * paramP); i++){
        int r = -1;
        while(used.count(r = randInt() % paramM)){}
        used.insert(r);
        list[i] = r;
    }
}

void getOneExchange(int *paths, int *ex){
    int *k = new int[paramM / (paramP * paramP)];
    randomExchangeList(k);
    for(int i = 0; i < paramM / (paramP * paramP); i++){
        for(int j = 0; j < paramN; j++){
            ex[i * paramN + j] = paths[k[i] * paramN + j];
        }
    }
    delete[] k;
}

void setOneExchange(int *paths, int *ex){
    unordered_set<int> used;
    int *k = new int[paramM / (paramP * paramP)];
    randomExchangeList(k);
    for(int i = 0; i < paramM / (paramP * paramP); i++){
        for(int j = 0; j < paramN; j++){
            paths[k[i] * paramN + j] = ex[i * paramN + j];
        }
    }
    delete[] k;
}

void getExchange(int *paths, int *ex){
    for(int i = 0; i < (paramM / (paramP * paramP)) * paramN; i += paramN){
        getOneExchange(paths, ex + i);
    }
}

void setExchange(int *paths, int *ex){
    for(int i = 0; i < (paramM / (paramP * paramP)) * paramN; i += paramN){
        setOneExchange(paths, ex + i);
    }
}

void doExchange(int *paths){
    int *exGet = new int[(paramM / (paramP * paramP)) * paramP * paramN];
    int *exSet = new int[(paramM / (paramP * paramP)) * paramP * paramN];
    int *len = new int[paramP];
    int *offset = new int[paramP];
    for(int i = 0; i < paramP; i++){
        len[i] = (paramM / (paramP * paramP)) * paramN;
        offset[i] = (paramM / (paramP * paramP)) * paramN * i;
    }

    getExchange(paths, exGet);
    MPI_Alltoallv(exGet,len,offset,MPI_INT,exSet,len,offset,MPI_INT,MPI_COMM_WORLD);
    getExchange(paths, exSet);
    delete[] exGet;
    delete[] exSet;
    delete[] len;
    delete[] offset;
}

int getGlobalMin(int min, int *minPath, int *gMinPath){
    int *sendMin = new int[paramP * (paramN + 1)];
    int *recvMin = new int[paramP * (paramN + 1)];
    int *len = new int[paramP];
    int *offset = new int[paramP];

    for(int i = 0; i < paramP; i++){
        sendMin[i * (paramN + 1)] = min;
        len[i] = paramN + 1;
        offset[i] = i * (paramN + 1);
        for(int j = 0; j < paramN; j++){
            sendMin[i * (paramN + 1) + j + 1] = minPath[j];
        }
    }

    //cout << "proc" << myID << " mins: " << min << endl;
    MPI_Alltoallv(sendMin,len,offset,MPI_INT,recvMin,len,offset,MPI_INT,MPI_COMM_WORLD);
    int gMin = recvMin[0];
    if(myID==0){cout << "current round mins: ";}
    for(int i = 0; i < paramP; i++){
        if(recvMin[i * (paramN + 1)] < gMin){
            gMin = recvMin[i * (paramN + 1)];
            for(int j = 0; j < paramN; j++){
                gMinPath[j] = recvMin[i * (paramN + 1) + j + 1];
            }
        }
        if(myID==0){cout << recvMin[i * (paramN + 1)] << " ";}
    }
    if(myID==0){cout << "with global min: " << gMin << endl;}

    delete[] sendMin;
    delete[] recvMin;
    delete[] len;
    delete[] offset;

    return gMin;
}

int main(int argc,char* argv[]){

    int id;
	int p;
	double wtime;
    
    initEverything(stoi(argv[1]), stoi(argv[2]), stoi(argv[3]), stoi(argv[4]), stof(argv[5]), stof(argv[6]), stof(argv[7]), stoi(argv[8]));
    //cout << paramN << " " << paramM << " " << paramK << " " << paramR << " " << paramS << " " << paramU << " " << paramZ << endl;

	MPI::Init(argc, argv); //  Initialize MPI.
	p = MPI::COMM_WORLD.Get_size(); //  Get the number of processes.
	myID = id = MPI::COMM_WORLD.Get_rank(); //  Get the individual process ID.

    //updateM(paramM/p);

    initRand(id);
    
    paramD = new int[paramN * paramN];
    paramP = p;
    
    if(id == 0){
        string fn(argv[9]);
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

    int gMin = -1;

    int * minPath = new int[paramN];
    int * gMinPath = new int[paramN];

    int countR = 0;
    int countK = 0;
    int countE = 0;

    while(true){
        cout << myID << "dist1" << endl;
        for(int i = 0; i < paramM; i++){
            distances[i] = roundDistance(paths + i * paramN);
        }
        cout << myID << "dist2" << endl;
        int currminindex = minDistanceIndex(paths)
        int currmin = distances[currminindex];

        cout << myID << "min1" << currmin << endl;
        if(min == -1 || currmin < min){
            min = currmin;
            for(int i = 0; i < paramN; i++){
                minPath[i] = paths[currminindex * paramN + i];
            }
        }
        cout << myID << "min2" << endl;
        
        countE = (countE + 1) % paramE;

        if(countE == 0){
            countK++;
            int t = getGlobalMin(min, minPath, gMinPath);
            if(t == gMin){
                countR++;
            }else{
                gMin = t;
                countR = 0;
            }
            if(countR >= paramR || countK >= paramK){
                break;
            }
            doExchange(paths);
        }

        //cout << "Current round min distance: " << currmin << ", all time min: " << min << ", keep for " << countR << " rounds." << endl;

        
        cout << myID << "sel1" << endl;
        selection(distances, survivors);
        cout << myID << "sel2" << endl;

        int *tempPtr = paths;
        paths = paths2;
        paths2 = tempPtr;

        cout << myID << "cp1" << endl;
        copySelection(paths2, paths, survivors);
        cout << myID << "cp2" << endl;

        cout << myID << "co1" << endl;
        fillCrossover(paths);
        cout << myID << "co2" << endl;
        cout << myID << "mut1" << endl;
        doMutation(paths);
        cout << myID << "mut2" << endl;
    }
    
    if(id == 0){
        cout << "min path: ";
        for(int i = 0; i < paramN; i++){
            cout << gMinPath[i] << " ";
        }
        cout << endl;
    }

    MPI::Finalize();
}
