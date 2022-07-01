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
           "-o <string> specify output file name    \n"
           "-v          verbose output              \n"
    );
}

int main(int argc, char** argv){
    string inFilename;
    string outFilename;

    // default parameters
    corrector.numIter  = 10;
    corrector.maxLevel = 3;
    corrector.verbose  = false;

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
                    corrector.numIter = n;
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
        else if( strcmp(argv[i], "-v"   ) == 0 ){
            corrector.verbose = true;
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
        printf("loaded %s\n", inFilename.c_str());
    }
    else{
        printf("error loading %s\n", inFilename.c_str());
        return -1;
    }

    corrector.Correct();

    if(!outFilename.empty()){
        loader.Save(outFilename, &pg);
        printf("output written to %s\n", outFilename.c_str());
    }

    return 0;
}
