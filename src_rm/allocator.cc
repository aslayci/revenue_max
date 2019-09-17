#include "allocator.h"
#include "anyoption.h"
#include "memoryusage.h"

namespace _Cide{
    
    allocator::allocator(AnyOption* opt1) {
        
        opt = opt1;
        delim = " \t";
        
        n = strToInt(opt->getValue("n"));
        m = strToInt(opt->getValue("m"));
        nrTopics = strToInt(opt->getValue("nrTopics"));
        nrCompanies = strToInt(opt->getValue("nrCompanies"));
        epsilon = strToFloat(opt->getValue("epsilon"));
        ell = strToFloat(opt->getValue("ell"));
        alpha = strToFloat(opt->getValue("alpha"));
        outFolderName = opt->getValue("outputFolder");
        windowSize = strToInt(opt->getValue("windowSize"));
        
        if(string(opt->getValue("greedyCriteria")).compare("ca") == 0) {
            windowSize = 1;
        }
        
        delta = pow(float (1 / (float) n), ell);
        
        for(int i = 0; i < n; i++) // keeps p^z_{uv} for each topic z
            graphT.push_back(vector<int>());
        
        advList = new advertiserList();
        //risList = new risEstimatorList();
        timList = new TimGraphList();
        
        // create advertiser objects and initialize their TIC prob-vector probT, aligned with graphT
        for(int i = 0; i < nrCompanies; i++) {
            advertiser *aa = new advertiser(i,nrTopics);
            advList->push_back(aa);
            _Cide::TimGraph *tim = new _Cide::TimGraph(aa, epsilon, n, m, string(opt->getValue("greedyCriteria")), windowSize);
            timList->push_back(tim);
            aa->currentPayment = 0;
            aa->currentRevenue = 0;
            aa->currentSeedCosts = 0;
            aa->maxCost = 0;
            aa->minCost = (float) n ;
            aa->seedUserCosts.resize(n);
            for (int i = 0; i < n; i++)
                aa->probT.push_back(std::vector< float>());
        }
        
        // advertiser object should be initialized before these calls
        readItemDistsFile();
        readTICGraph();
        attentionQuota = std::vector<int>(n,1); // only disjointness
        readBudgetsFile();
        readIncentiveCosts(); // reading incentive costs as matrix nodes X ads
        arrangeOutputFiles();
        allocateSeeds();
        cout << "finished allocating seeds.." << endl;
        
    }
    
    // allocator for revmax
    void allocator::allocateSeeds() {
        
        allocQueue.clear();
        
        advertiser *adv;
        TimGraph *tim;
        infPair bests;
        
        cout << string(opt->getValue("greedyCriteria")) << " allocation started: using windows =  " << windowSize << endl;
        time(&startTime);	// get it for the total algo running time
        
        // initial RR sets computation and allocQueue update
        for(int iterC = 0; iterC < nrCompanies; iterC++) {
            adv = advList->at(iterC);
            tim = timList->at(iterC);
//            cout << "adv max cost: " << adv->maxCost << endl;
            tim->doInitialGeneration();
            
            // cost-agnostic
            if(string(opt->getValue("greedyCriteria")).compare("ca") == 0)
                tim->findBestCANode(bests);
            
            // cost-sensitive
            else
                tim->findBestCSNode(bests);
            
            allocQueue.insert(pair< float,int>(bests.second, iterC));
            
        }
        
        
        // iterated allocation starts
        while(!allocQueue.empty()) {	//stop allocation when no more advertiser is available for allocation
            
            // give allocation priority to the advertiser who has the highest value of criter
            multimap< float, int>::iterator i = allocQueue.end();
            i--;
            
            adv = advList->at(i->second);
            tim = timList->at(i->second);
            
            allocQueue.erase(i); //remove from allocation queue, next best will be reinserted soon if allocation do not finish for this
            
            if(attentionQuota[tim->candidateNode] > 0) { // if the best candidate still has attention quota left for assignment
                
                if((adv->currentPayment + tim->candidateMgRev + (adv->seedUserCosts[tim->candidateNode])) <= adv->budget) {
                    
                    tim->assignBestNode(); // inserts into advertiser's seed set and decreases the node's attention quota
                    adv->currentPayment += (tim->candidateMgRev + adv->seedUserCosts[tim->candidateNode]);
                    adv->currentRevenue += tim->candidateMgRev;
                    adv->currentSeedCosts += adv->seedUserCosts[tim->candidateNode];
                    
                    if((int) tim->seedSet.size() == tim->kappa) {
                        tim->updateEstimates();
                        //cout << "updating estimates for advertiser " << adv->advertiserId << endl;
                    }
                    
                    if(string(opt->getValue("greedyCriteria")).compare("ca") == 0) // cost-agnostic
                        tim->findBestCANode(bests);
                    
                    else // cost-sensitive
                        tim->findBestCSNode(bests);
                    
                    allocQueue.insert(pair< float,int>(bests.second, adv->advertiserID));
                    continue;
                    
                } // end of budget feasible part
                
                else{ // if not budget feasible, we stop the allocation for this advertiser so sktir et
                    //cout << "budget exhausted " << endl;
                    continue; 
                }
                
            }
            
            else { // select another best node for this advertiser cause this seed is already allocated to another advertiser due to allocation priority
                
                if(string(opt->getValue("greedyCriteria")).compare("ca") == 0)
                    tim->findBestCANode(bests);
                
                else
                    tim->findBestCSNode(bests);
                
                allocQueue.insert(pair< float,int>(bests.second, adv->advertiserID));
                
                continue;
                
            }
            
            
        } // while in sonu, alloc queue empty buradan sonra
        
        
        float totalDuration = getRunningTime(startTime); // in seconds
        float totalMemory = disp_mem_usage(); // in MB
        cout << "total time taken (in seconds) " << totalDuration << "total memory (in mb)" << totalMemory << endl;
        
        
        // write results to master and adv specific files
        float total_revenue = 0;
        float total_payment = 0;
        float total_budget = 0;
        int seedSizeTotal = 0;
        float totalIncentiveCosts = 0;
        
        
        cout << "writing the results to output files.." << endl;
        
        for(int i = 0; i < nrCompanies; i++) {
            adv = advList->at(i);
            tim = timList->at(i);
            
            for (int j = 0; j < (int) tim->seedSet.size(); j++) {
                writeInAdvertisersFile(i, tim->seedSet[j], tim->seedMgRevs[j],(adv->seedUserCosts[tim->seedSet[j]]));
            }
            
            cout << "kontrol : payment and revenue " << adv->currentPayment << " & " << adv->currentRevenue << endl;
            
            writeInMasterOutputFile(algorithmTip, costFuncTip, windowTip, alpha, windowSize, adv->advertiserID, adv->currentRevenue, adv->currentSeedCosts, adv->currentPayment, adv->budget, tim->seedSet.size(), totalDuration, totalMemory);
            
            total_revenue += adv->currentRevenue;
            total_payment += adv->currentPayment;
            seedSizeTotal += tim->seedSet.size();
            totalIncentiveCosts += adv->currentSeedCosts;
            total_budget += adv->budget;
        }
        
        
        writeInMasterOutputFile(algorithmTip, costFuncTip, windowTip, alpha, windowSize, -1, total_revenue, totalIncentiveCosts, total_payment, total_budget, seedSizeTotal, totalDuration, totalMemory);
        
        cout << "master file header : " << endl;
        cout << "algorithm costFunctionType windowTip alpha window advID revenue seedCosts payment budget seedSetSize runTime(sec) memory(mb)" << endl;
        
        
        
    }
    
    
    void allocator::readTICGraph() {
        
        string probGraphFile = opt->getValue("probGraphFile");
        cout << "Reading file " << probGraphFile << endl;
        ifstream myfile (probGraphFile.c_str(), ios::in);
        
        float *dists;
        float p;
        advertiser *advTemp;
        
        int nrEdges = 0;
        set<int> nodes; // kontrol amacli simdi
        
        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;
                nrEdges++;
                
                std::string::size_type pos = line.find_first_of(delim);
                int prevpos = 0;
                
                //first user
                string str = line.substr(prevpos, pos-prevpos);
                int u1 = strToInt(str);
                
                //second user
                prevpos = line.find_first_not_of(delim, pos);
                pos = line.find_first_of(delim, prevpos);
                int u2 = strToInt(line.substr(prevpos, pos-prevpos));
                
                if (u1 == u2)
                    continue;
                
                graphT[u2].push_back(u1); //insert to the transposed graph
                
                // kontrol amacli
                nodes.insert(u1);
                nodes.insert(u2);
                
                prevpos = line.find_first_not_of(delim, pos);
                
                str = line.substr(prevpos);
                dists = new float[nrTopics];
                stringTokenizer(str, dists, nrTopics, delim);
                
                for(int i = 0; i < nrCompanies; i++) {
                    advTemp = advList->at(i);
                    p = 0.0;
                    for(int j = 0; j < nrTopics; j++)
                        p += (dists[j] * advTemp->gamma[j]);
                    
                    advTemp->probT[u2].push_back(p);
                }
            }
            
            myfile.close();
        }
        
        else
            cout << "Can't open friendship graph file " << probGraphFile << endl;
        
        cout << "Built transposed adj matrix from file " << endl;
        cout << "number of nodes " << nodes.size() << endl;
        cout << "number of edges " << nrEdges << endl;
        
    }
    
    void allocator::readItemDistsFile() {
        
        cout << "reading item distributions file " << endl;
        string itemDistsFile = opt->getValue("itemDistsFile");
        ifstream myfile(itemDistsFile.c_str(), ios::in);
        float *tokens;
        
        int advIndex = 0;
        
        if(myfile.is_open()) {
            while(!myfile.eof()) {
                std::string line;
                getline(myfile, line);
                if(line.empty())
                    continue;
                tokens = new  float[nrTopics];
                stringTokenizer(line, tokens, nrTopics, delim);
                advList->at(advIndex++)->setItemDist(tokens, nrTopics);
                if(advIndex >= nrCompanies)
                    break;
            }
            
            myfile.close();
        }
        
        else {
            cout << "problem opening the item distributions file, exiting... " << itemDistsFile <<  endl;
            exit(1);
        }
        
    }
    
    // this is just to create individual ic graph files from tic
    //    void allocator::writeICGraphFile() {
    //
    //        string fileName_ = "input/flixster/Flixster_ICGraph_item_";
    //        string extension = ".txt";
    //        advertiser *adv;
    //
    //        for (int i = 0; i < nrCompanies; i++) {
    //            string fileName = fileName_ + intToStr(i) + extension;
    //            cout << "working on file " << fileName << endl;
    //            ofstream outStream;
    //            if(outStream.is_open())
    //                outStream.close();
    //            outStream.open(fileName.c_str());
    //            adv = advList->at(i);
    //
    //            for (int j = 0; j < n; j++) {
    //                for (int hede = 0; hede < (int) graphT[j].size(); hede++) {
    //                    outStream << graphT[j][hede] << "\t" << j << "\t" << adv->probT[j][hede] << endl;
    //                }
    //            }
    //            outStream.close();
    //        }
    //
    //    }
    
    
    void allocator::readBudgetsFile() {
        
        cout << "reading advertisers' budgets and cpes " << endl;
        string budgetsFile = opt->getValue("budgetsFile");
        ifstream myfile(budgetsFile.c_str(), ios::in);
        float *bc;
        
        advertiser *adv;
        int advIndex = 0;
        if(myfile.is_open()) {
            while(!myfile.eof()) {
                std::string line;
                
                getline(myfile, line);
                if(line.empty())
                    continue;
                
                bc = new  float[2];
                stringTokenizer(line, bc, 2, delim);
                
                adv = advList->at(advIndex++);
                adv->budget = bc[0];
                adv->cpe = bc[1];
                
                if(advIndex >= nrCompanies)
                    break;
            }
            
            myfile.close();
        }
        
        else {
            cout << "problem opening the budgets file, exiting... " << budgetsFile <<  endl;
            exit(1);
        }
        
        // kontrol
        cout << "budget kontrol " << endl;
        for(int advIndex = 0; advIndex < nrCompanies; advIndex++)
            cout << advList->at(advIndex)->budget << endl;
        
        cout << "cpe kontrol " << endl;
        for(int advIndex = 0; advIndex < nrCompanies; advIndex++)
            cout << advList->at(advIndex)->cpe << endl;
    }
    
    // this reads a cost file in the form of a cost matrix in the form of nodes X ads -- might be different for the scalability version
    void allocator::readIncentiveCosts() {
        
        string readIncentiveCostsFile = opt->getValue("incentiveCostsFile");
        ifstream myfile(readIncentiveCostsFile.c_str(), ios::in);
        
        int lineIndex = 0;
        float *tokens;
        float tempCostToken = 0.0;
        
        if(myfile.is_open()) {
            while(!myfile.eof()) {
                std::string line;
                getline(myfile, line);
                if(line.empty())
                    continue;
                tokens = new float[nrCompanies];
                stringTokenizer(line, tokens, nrCompanies, delim);
                
                for (int i = 0; i < nrCompanies; i++) {
                    if(tokens[i] < 1.0)
                        tokens[i] = 1.1;
                }
                
                for(int i = 0; i < nrCompanies; i++) {
                    
                    if(string(opt->getValue("costFunctionType")).compare("l") == 0) { // linear
                        tempCostToken = tokens[i] * alpha;
                    }
                    
                    else if(string(opt->getValue("costFunctionType")).compare("u") == 0) { // uniform -- reads uniform input
                        tempCostToken = tokens[i] * alpha;
                    }

                    else if(string(opt->getValue("costFunctionType")).compare("q") == 0) { // quadratic
                        tempCostToken = tokens[i] * tokens[i] * alpha;
                    }
                    
                    else if(string(opt->getValue("costFunctionType")).compare("s") == 0) { // sublinear
                        tempCostToken = log(tokens[i]) * alpha;
                    }
                    
                    else if(string(opt->getValue("costFunctionType")).compare("r") == 0) { // random
                        
                        tempCostToken = tokens[i] * alpha; // daha sonra shuffle edilecek bu
                    }
                    
                    advList->at(i)->seedUserCosts[lineIndex] = tempCostToken;
                    if(advList->at(i)->maxCost < tempCostToken)
                        advList->at(i)->maxCost = tempCostToken;
                    
                    if(advList->at(i)->minCost > tempCostToken)
                        advList->at(i)->minCost = tempCostToken;                    
                }
                
                lineIndex++;
                
            }
            
            myfile.close();
            
            //            cout << " coost - kontrol : " << advList->at(0)->seedUserCosts[0] << endl;
            
            // for random cost function shuffle the vectors here
            if(string(opt->getValue("costFunctionType")).compare("r") == 0) { // random
                for (int i =0 ; i < nrCompanies; i++) {
                    //cout << "adv " << i << " cost kontrol " << advList->at(i)->maxCost << endl;
                    //                    cout << "randomness kontrol b4 : " << advList->at(i)->seedUserCosts[0] << endl;
			std::srand(std::time(0));
                    std::random_shuffle(advList->at(i)->seedUserCosts.begin(), advList->at(i)->seedUserCosts.end());
                    //                    cout << "randomness kontrol after : " << advList->at(i)->seedUserCosts[0] << endl;
                }
            }
            
            
            cout << "cost control " << endl;
            for (int i = 0; i < nrCompanies; i++) {
                cout << "max for adv " << i << " " << advList->at(i)->maxCost << endl;
                cout << "min for adv " << i << " " << advList->at(i)->minCost << endl;
            }
            
            
        }
        
        else {
            cout << "problem opening the incentive costs file, exiting... " << readIncentiveCostsFile <<  endl;
            exit(1);
        }
        
    }
    
    allocator::~allocator() {
        //cout << "destructor called " << endl;
    }
    
    void allocator::arrangeOutputFiles() {
        
        string command = string("mkdir -p ") + outFolderName ;
        system(command.c_str());
        string fname;
        
        if(string(opt->getValue("costFunctionType")).compare("l") == 0)
            costFuncTip = "linear";
        
        else if(string(opt->getValue("costFunctionType")).compare("u") == 0)
            costFuncTip = "uniform";

        else if(string(opt->getValue("costFunctionType")).compare("q") == 0)
            costFuncTip = "quadratic";
        
        else if(string(opt->getValue("costFunctionType")).compare("s") == 0)
            costFuncTip = "sublinear";
        
        else if(string(opt->getValue("costFunctionType")).compare("r") == 0)
            costFuncTip = "random";
        
//        costFuncTip = costFuncTip + floatToStr(alpha);
        string costFuncTip2 = costFuncTip + floatToStr(alpha);
        
        if(string(opt->getValue("greedyCriteria")).compare("ca") == 0) {
            algorithmTip = "TI-CARM";
            windowTip = "NA";
            fname = "ca_" + costFuncTip2 + "_w" + windowTip;
        }
        
        else if(string(opt->getValue("greedyCriteria")).compare("cs") == 0) {
            algorithmTip = "TI-CSRM";
            if(windowSize == 0) windowTip = "NA";
            else windowTip = intToStr(windowSize);
            fname = "cs_" + costFuncTip2 + "_w" + windowTip;
            
        }
        
        string masterFileName = "tirm_master_" + fname + ".txt";

        //cout << masterFileName << endl;
        
        /* arrange the output files - folders */
        outMasterName = outFolderName + OS_SEP + masterFileName;
        openOutputFiles();
        //cout << "opened master " << endl;
        
        //
        outFileStreams = new ofstream[nrCompanies];
        outFileNames = new string[nrCompanies];
        /* create output file names for each advertiser */
        for(int i =0; i < nrCompanies; i++) {
            outFileNames[i] = outFolderName + OS_SEP + "tirm_adv_" + intToStr(i) + "_" + fname  + ".txt";
            cout << outFileNames[i] << endl;
        }
        
        for(int i =0; i < nrCompanies; i++) {
            openOutputFiles(advList->at(i)->advertiserID);
        }
        //
        // memory icin
        command = string("mkdir -p temp") ;
        system(command.c_str());
    }
    
    void allocator::writeInAdvertisersFile(int advID, int v, float mgRev, float seedCost) {
        (outFileStreams[advID]) << v << " " << mgRev << " " << seedCost << endl;
    }
    
    void allocator::openOutputFiles(int advID) {
        
        // will create output files per advertiser
        
        if((outFileStreams[advID]).is_open())
            (outFileStreams[advID]).close();
        
        (outFileStreams[advID]).open(outFileNames[advID].c_str());
        
        if ((outFileStreams[advID]).is_open() == false) {
            cout << "Can't open file " << outFileNames[advID]  << " for writing" << endl;
            exit(1);
        }
    }
    
    void allocator::openOutputFiles() {
        // will create master-output file
        if(outMasterStream.is_open())
            outMasterStream.close();
        
        outMasterStream.open(outMasterName.c_str());
        
        if (outMasterStream.is_open() == false) {
            cout << "Can't open file " << outMasterName  << " for writing" << endl;
            exit(1);
        }
        
//        outMasterStream << "algorithm costFunctionType windowTip alpha window advID revenue seedCosts payment budget seedSetSize runTime(sec) memory(mb)" << endl;
        
    }
    
    
    void allocator::writeInMasterOutputFile(string algorithmTip, string costFuncTip, string windowTip, float alpha, int windowSize, int advertiserID, float revenue, float seedCosts, float payment, float budget, int size, float duration, float memory) {
        
        // algorithm costFunctionType windowTip alpha window advID revenue seedCosts payment budget seedSetSize runTime(sec) memory(mb)
        
        outMasterStream << algorithmTip << " " << costFuncTip << " " <<  windowTip << " " << alpha << " " << windowSize << " " << advertiserID << " " << revenue << " " << seedCosts << " " << payment << " " << budget << " " << size << " " << duration << " " << memory << endl;
    }

    
    
    
    
}
