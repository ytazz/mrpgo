#include <mrpgo.h>
#include <windows.h>

using namespace std;
using namespace mrpgo;

Posegraph pg;
Loader    loader;
Corrector corrector(&pg);

void PrintHelp(){
    printf("Usage: mrpgo [options] [filename]       \n"
           "                                        \n"
           "General options:                        \n"
           "----------------------------------------\n"
           "-help / -h           Displays this help.\n"
           "                                        \n"
           "Program Options:                        \n"
           "----------------------------------------\n"
           "-i <int>    specify number of iterations\n"
           "-L <int>    specify max level           \n"
           "-s <string> specify statistics file name\n"
           "-o <string> specify output file name    \n"
           "-v          verbose output              \n"
    );
}

int main(int argc, char** argv){
    string inFilename;
    string outFilename;

    // default parameters
    int numIter = 10;
    corrector.maxLevel = 3;
    //corrector.solverType     = Corrector::SolverType::Cholmod;
    corrector.solverType     = Corrector::SolverType::Custom;
    corrector.param.verbose  = false;

    // process command line arguments
    if(argc == 1){
        PrintHelp();
        return 0;
    }
    for(int i = 1; i < argc; i++){
        if( strcmp(argv[i], "-h"   ) == 0 || 
            strcmp(argv[i], "-help") == 0 ){
            PrintHelp();
        }
        else if( strcmp(argv[i], "-i"   ) == 0 ){
            if(++i < argc){
                int n = atoi(argv[i]);
                if(n > 0){
                    numIter = n;
                }
                else{
                    printf("error: number of iterations must be > 0\n");
                }
            }
        }
        else if( strcmp(argv[i], "-L"   ) == 0 ){
            if(++i < argc){
                int L = atoi(argv[i]);
                if(L >= 0){
                    corrector.maxLevel = L;
                }
                else{
                    printf("error: number of iterations must be >= 0\n");
                }
            }
        }
        else if( strcmp(argv[i], "-p"   ) == 0 ){
            if(++i < argc){
                corrector.numThreads = atoi(argv[i]);
            }
        }
        else if( strcmp(argv[i], "-v"   ) == 0 ){
            corrector.param.verbose = true;
        }
        else if( strcmp(argv[i], "-s"   ) == 0 ){
            if(++i < argc){
                corrector.statFilename = argv[i];
            }
        }
        else if( strcmp(argv[i], "-o"   ) == 0 ){
            if(++i < argc){
                outFilename = argv[i];
            }
        }
        else{
            inFilename = argv[i];
        }
    }

    if(inFilename.empty()){
        printf("filename not specified\n");
        return -1;
    }

    if(loader.Load(inFilename, &pg)){
        printf("loaded %s as %s posegraph\n", inFilename.c_str(), pg.space == Posegraph::Space::SE2 ? "2D" : "3D");
    }
    else{
        printf("error loading %s\n", inFilename.c_str());
        return -1;
    }

    for(int i = 0; i < numIter; i++){
        corrector.Correct();
    }

    if(!outFilename.empty()){
        loader.Save(outFilename, &pg);
        printf("output written to %s\n", outFilename.c_str());
    }

    return 0;
}
