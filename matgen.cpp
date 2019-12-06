#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <math.h>
#include <time.h>

using namespace std;

int paramN = 100; //# of sites
int maxDist = 10000;

int main(int argc,char* argv[]){
    string fn(argv[1]);
    ofstream out(fn);
    srand(time(NULL));

    paramN = stoi(argv[2]);
    maxDist = stoi(argv[3]);

    int *mat = new int[paramN * paramN];

    for(int i = 0; i < paramN; i++){
        for(int j = 0; j < paramN; j++){
            if(i == j){
                mat[i * paramN + j] = 0;
            }else if(j > i){
                mat[i * paramN + j] = rand() % maxDist;
            }else{
                mat[i * paramN + j] = mat[j * paramN + i];
            }
        }
    }

    for(int i = 0; i < paramN; i++){
        for(int j = 0; j < paramN; j++){
            out << mat[i * paramN + j] << " ";
        }
        out << endl;
    }
    out.close();
}