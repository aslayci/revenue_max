#ifndef ALLOCATOR_H
#define ALLOCATOR_H


#include <ctime>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include "anyoption.h"
#include "utils.h"
#include "advertiser.h"
#include "TimGraph.h"
#include <cmath>
#include <algorithm>



namespace _Cide{
    
    typedef std::vector<advertiser*> advertiserList;
    //typedef std::vector<risEstimator*> risEstimatorList;
    typedef std::vector<TimGraph*> TimGraphList;
    
    class allocator {
        
    public:
        
        time_t startTime;
    
        AnyOption *opt;
        int n, m, nrTopics, nrCompanies;
        float epsilon, delta, alpha, ell;
        string delim, probGraphFile, outputCostsFile;
        
        advertiserList *advList;
        //risEstimatorList *risList;
        TimGraphList *timList;
        
        allocator(AnyOption* opt);
        ~allocator();
        
        multimap< float,int> allocQueue; // allocation priority
        
        
        string algorithmTip; // for plots
        string costFuncTip; // for plots
        string windowTip; // for plots (string cause will define "NA")
        int windowSize; 
        
        
        
        void readTICGraph();
        void readItemDistsFile();
        void readBudgetsFile();
        void readSeedSetForSimulation(std::vector<int> &seedSetForSim); // bunu verify etmek icin kullanabilirim
        void readIncentiveCosts();
        void allocateSeeds();
        
        void arrangeOutputFiles();
        void openOutputFiles(int advID);
        void openOutputFiles();
        void writeInAdvertisersFile(int advID, int v, float mgRev, float seedCost);
        void writeInMasterOutputFile(string algorithmTip, string costFuncTip, string windowTip, float alpha, int windowSize, int advertiserID, float revenue, float seedCosts, float payment, float budget, int size, float duration, float memory);
            
        
        ofstream *outFileStreams;
        ofstream outMasterStream;
        string outFolderName;
        string *outFileNames; // array of output file names for each advertiser
        string outMasterName;
        
        void writeICGraphFile();
        
    
    
    
    };


}


#endif
