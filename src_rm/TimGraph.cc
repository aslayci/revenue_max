#ifndef _TIMGRAPH_C_
#define _TIMGRAPH_C_

#include "TimGraph.h"

namespace _Cide{
        
    TimGraph::TimGraph(advertiser *adv,  float eps, int nrNodes, int NrEdges, string greedyKriter, int windowSize) {
        
        // 		cout << "........ Initializing TimGraph object ............" << endl;
        // initialize vars
        this->epsilon = eps;
        this->adv = adv;
        this->n = nrNodes;
        this->m = NrEdges;
        this->seedSet.clear();
        
        
        this->windowSize = windowSize;
        
        if(greedyKriter.compare("ca") == 0) {
            this->criterGreedy = criterGreedy;
            windowSize = 1; // ca this is by default
        }
        
        // branch the cost-sensitive between no windows and windows
        else if(greedyKriter.compare("cs") == 0 && windowSize == 0) { // no windows
            this->criterGreedy = "csn"; // normal cost-sensitive
        }
        
        else { // cost-sensitive with windows
            this->criterGreedy = "csw";
        }
        
        
        visit_mark = std::vector<int>(n,0);
        visit = std::vector<bool>(n,false);
        hyper_degree = std::vector<int>(n,0);
        //nodeCTRs = std::vector<float>(n,0);
        
        for(int i = 0; i < 12; i++)
            sfmt_init_gen_rand(&sfmtSeed, i+1234);
        
    }
    
    TimGraph::~TimGraph(void) {
        cout << "destructor called for timgraph" << endl;
    }
    
    void TimGraph::doInitialGeneration() {
        
        usersExamined.clear();
        criterQueue.clear();
        
        kappa = 1;
        
        theta = theta_old = 0;
        estimateTheta();
        cout << "advertiser " << adv->advertiserID << " estimated theta " << theta << endl;
        
        for(int i = 0; i < theta; i++) {
            hyperGT_adv.push_back(std::vector<int>());
            isCovered.push_back(false);
            BuildHypergraphNode(sfmt_genrand_uint32(&sfmtSeed)%n, i, true, hyperGT_adv);
        }
        
        hyperG_adv.clear();
        for(int i = 0; i < n; i++)
            hyperG_adv.push_back(std::vector<int>());
        
        for(int i = 0; i < theta; i++){
            for(int j = 0; j < (int) hyperGT_adv[i].size(); j++) {
                int t = hyperGT_adv[i][j];
                hyperG_adv[t].push_back(i);
            }
        }
        
        
        for(int i = 0; i < n; i++) {
            
            hyper_degree[i] = hyperG_adv[i].size();
            
            if(criterGreedy.compare("csn") == 0) { // cost-sensitive with no windows
                mmg_temp = (( float) n * (( float) hyper_degree[i] / theta)) * adv->cpe;
                criterQueue.insert(pair<float,int>(mmg_temp / adv->seedUserCosts[i],i));
            }
            
        }
        
    }
    
    // returns the best cost-agnostic node -- without optimization
    void TimGraph::findBestCANode(infPair &bestCandidate) {
        
        float mmg;
        
        int id = -1;
        float maxVal = 0.0;
        
        for(int i = 0; i < (int) hyper_degree.size(); i++) { // here "i" is the node id
            if((attentionQuota[i] > 0) && (maxVal < hyper_degree[i])) {
                maxVal = hyper_degree[i];
                id = i;
            }
        }
        
        mmg = (( float) n * (( float) maxVal / theta)) * adv->cpe;
        this->candidateNode = id;
        this->candidateNrRR = hyper_degree[id];
        this->candidateMgRev = mmg;
        bestCandidate = std::make_pair(id,mmg);
        hyper_degree[id]=0;  // assign zero since either it will be allocated soon or will be out of the game due to attention bound

        
    }
    
    
    // returns the best cost-sensitive node -- with window
    void TimGraph::findBestCSNode(infPair &bestCandidate) {
        
        if(windowSize == 0) { // no windows
            
            multimap<float,int>::iterator it = criterQueue.end();
            it--;
            
            while (true) {
                
                if(attentionQuota[it->second] == 0) {
                    criterQueue.erase(it);
                    it = criterQueue.end();
                    it--;
                    continue;
                }
                
                if(usersExamined.find(it->second) == usersExamined.end()) { //this node is not explored after assignBestNode (ie. hyper_degree unchanged)
                    
                    this->candidateNode = it->second;
                    this->candidateNrRR = hyper_degree[it->second];
                    this->candidateMgRev = it->first * adv->seedUserCosts[this->candidateNode];
                    bestCandidate = std::make_pair(this->candidateNode,it->first); // should contain the ratio for CS
                    hyper_degree[candidateNode]=0;  // assign zero since either it will be allocated soon or will be out of the game due to
                    criterQueue.erase(it);
                    break;
                    return;
                    
                }
                
                
                else {
                    int idTemp = it->second;
                    mmg_temp = (( float) n * (( float) hyper_degree[idTemp] / (float) theta)) * adv->cpe;
                    criterQueue.erase(it);
                    usersExamined.erase(idTemp);
                    criterQueue.insert(pair<float,int>(mmg_temp / adv->seedUserCosts[idTemp] ,idTemp)); // for CS, we put the ratio
                    it = criterQueue.end();
                    it--;
                }
            }
            
        }
        
        else { // with windows
        
            // cost-agnostic gibi ama top yerine top-windowSize tut
            
            float *windowMmgs = new float[windowSize];
            int *windowNodeIDs = new int[windowSize];
            float mmg = 0;
            int iter = 0;
            
            for (int i = 0; i < windowSize; i++) {
                windowMmgs[i] = 0.0;
                windowNodeIDs[i] = -1;
            }
            
            for(int i = 0; i < (int) hyper_degree.size(); i++) { // here "i" is the node id
                
                if(attentionQuota[i] == 0)
                    continue;
                
                mmg = (( float) n * (( float) hyper_degree[i] / theta)) * adv->cpe;
                
                if(windowNodeIDs[windowSize - 1] == -1 || mmg > windowMmgs[windowSize - 1]) {
                    
                    iter = windowSize - 1;
                    while(iter > 0 && mmg > windowMmgs[iter - 1]) {
                        windowNodeIDs[iter] = windowNodeIDs[iter - 1];
                        windowMmgs[iter] = windowMmgs[iter - 1];
                        iter--;
                    }
                    
                    if(mmg > windowMmgs[iter]) {
                        windowNodeIDs[iter] = i;
                        windowMmgs[iter] = mmg;
                    }
                
                }

            }
            
            // now check for the best cost-sensitive node among the nodes that reside in windowNodeIDs
            float maxCSVal = 0;
            float csKriter = 0;
            int id = -1;
            
            for (iter = 0; iter < windowSize; iter++) {
                csKriter = windowMmgs[iter] / adv->seedUserCosts[windowNodeIDs[iter]];
                if(maxCSVal < csKriter) {
                    maxCSVal = csKriter;
                    id = windowNodeIDs[iter];
                }
            }
            
            this->candidateNode = id;
            this->candidateNrRR = hyper_degree[id];
            this->candidateMgRev = (( float) n * (( float) hyper_degree[id] / theta)) * adv->cpe;
            bestCandidate = std::make_pair(id, maxCSVal);
            hyper_degree[id]=0;  // assign zero since either it will be allocated soon or will be out of the game due to attention
        
        }
            
    }
    
    //    // returns the best cost-sensitive node -- without window
    //    void TimGraph::findBestCSNode(infPair &bestCandidate) {
    //
    //        multimap<float,int>::iterator it = criterQueue.end();
    //        it--;
    //
    //        while (true) {
    //
    //            if(attentionQuota[it->second] == 0) {
    //                criterQueue.erase(it);
    //                it = criterQueue.end();
    //                it--;
    //                continue;
    //            }
    //
    //            if(usersExamined.find(it->second) == usersExamined.end()) { //this node is not explored after assignBestNode (ie. hyper_degree unchanged)
    //
    //                this->candidateNode = it->second;
    //                this->candidateNrRR = hyper_degree[it->second];
    //                this->candidateMgRev = it->first * adv->seedUserCosts[this->candidateNode]; // bu kisim cost-sensitive de carpi cost olarak gelecek
    //                bestCandidate = std::make_pair(this->candidateNode,this->candidateMgRev);
    //                hyper_degree[candidateNode]=0;  // assign zero since either it will be allocated soon or will be out of the game due to
    //                criterQueue.erase(it);
    //                break;
    //                return;
    //
    //            }
    //
    //
    //            else {
    //                int idTemp = it->second;
    //                mmg_temp = (( float) n * (( float) hyper_degree[idTemp] / (float) theta)) * adv->cpe;
    //                criterQueue.erase(it);
    //                usersExamined.erase(idTemp);
    //                criterQueue.insert(pair<float,int>(mmg_temp / adv->seedUserCosts[idTemp] ,idTemp)); // for CS, we put the ratio
    //                it = criterQueue.end();
    //                it--;
    //            }
    //        }
    //
    //    }
    
    void TimGraph::assignBestNode() {
        
        attentionQuota[candidateNode] = 0;
        hyper_degree[candidateNode] = 0;
        seedSet.push_back(candidateNode);
        nrRRCovered.push_back(candidateNrRR); // needed for spread estimate updates when the latent seed set size changes
        seedMgRevs.push_back(candidateMgRev);
        
        for(int j = 0; j < (int) hyperG_adv[candidateNode].size(); j++){ //for each RR set t covered by this node
            int t = hyperG_adv[candidateNode][j]; // rr set t
            if(!isCovered[t]){ //if the RR set was not covered before - since we do not delete covered from hyperG_adv except new RR set genration
                isCovered[t]=true; // bestnode covers it
                for(int z = 0; z < (int) hyperGT_adv[t].size(); z++){ //decrease the degree of parents who are associated with this RR set
                    int item = hyperGT_adv[t][z];
                    hyper_degree[item]--;
                    if(criterGreedy.compare("csn") == 0) // this is only for CS optimization
                        usersExamined.insert(item);
                    if(hyper_degree[item] < 0) { //myopic icin
                        hyper_degree[item] = 0;
                    }
                }
            }
        }
    }
    
    
    void TimGraph::updateEstimates() { // the seed nodes selected "by default" covers the new rr sets in their seed-selection order
        
        kappa += ceil((adv->budget - adv->currentPayment) / (candidateMgRev + adv->maxCost));
        
//        if(criterGreedy.compare("ca") == 0)
//            kappa += ceil((adv->budget - adv->currentPayment) / (candidateMgRev + adv->maxCost));
//        
//        // cide optim ski
//        else
//            kappa += ceil((adv->budget - adv->currentPayment) / (candidateMgRev + adv->minCost));
        
        //        cout << "adv: " << adv->advertiserID <<  " new latent seed set size " << kappa << endl;
        //        cout << "current reveneue wo update: " << adv->currentRevenue << endl;
        
        theta_old = theta;
        estimateTheta();
        
        if(theta <= theta_old) { // no updates, use the same set of RR sets
            theta = theta_old;
            //            cout << "theta remains the same " << endl;
        }
        
        else { // generate new additional RR sets and updates estimates
            
            //            cout << "new theta : " << theta << endl;
            
            // allocate more memory
            for(int i = theta_old; i < theta; i++) {
                hyperGT_adv.push_back(std::vector<int>());
                isCovered.push_back(false); // resize isCovered for new theta with false value
                BuildHypergraphNode(sfmt_genrand_uint32(&sfmtSeed)%n, i, true, hyperGT_adv);
            }
            
            
            // eskisi gibi yap
            hyperG_adv.clear();
            for(int i = 0; i < n; i++)
                hyperG_adv.push_back(std::vector<int>());
            
            for(int i = 0; i < theta; i++){
                if(!isCovered[i]) { //  we do not keep the RR sets that have been already covered before -- so that hyper_degree keeps the nr of uncovered RR sets
                    for(int j = 0; j < (int) hyperGT_adv[i].size(); j++) {
                        int t = hyperGT_adv[i][j];
                        hyperG_adv[t].push_back(i);
                    }
                }
            }
            
            for(int i = 0; i < n; i++)
                hyper_degree[i] = hyperG_adv[i].size();
            
            int oldNR;
            float totNR = 0; //ctr icin
            //                float costs = 0;
            
            for (int i = 0; i < (int) seedSet.size(); i++) {
                oldNR = nrRRCovered[i];
                nrRRCovered[i] = oldNR + hyper_degree[seedSet[i]];
                
                totNR += (float) nrRRCovered[i];
                //costs += (adv->seedUserCosts[seedSet[i]]);
                
                hyper_degree[seedSet[i]] = 0; // since we dont want to select this node as seed again
                
                // oldCPP
                for(int j = 0; j < (int) hyperG_adv[seedSet[i]].size(); j++){ //for each RR set covered by this node
                    int t = hyperG_adv[seedSet[i]][j]; // this is the RR set covered by this node
                    if(!isCovered[t]){
                        isCovered[t]=true; // thos old seed covers it
                        for(int z = 0; z < (int) hyperGT_adv[t].size(); z++){ //decrease the degree of parents who are associated with this RR set
                            int item = hyperGT_adv[t][z];
                            hyper_degree[item]--;
                            if(hyper_degree[item] < 0) {
                                hyper_degree[item] = 0;
                            }
                        }
                    }
                }
            }
            
            // update the revenue and the payment of advertisers wrt the new spread estimations
            adv->currentRevenue = (( float) n * (( float) totNR / theta)) * adv->cpe;
            adv->currentPayment = adv->currentRevenue + adv->currentSeedCosts;
            
            // update the local criteria queue -- only for cost-sensitive
            if(criterGreedy.compare("csn") == 0) {
                
                usersExamined.clear();
                criterQueue.clear();
                
                for(int i = 0; i < n; i++) {
                    
                    if(attentionQuota[i] == 0)
                        continue;
                    
                    mmg_temp = (( float) n * (( float) hyper_degree[i] / theta)) * adv->cpe;
                    criterQueue.insert(pair<float,int>(mmg_temp / adv->seedUserCosts[i],i));
                    
                }
            }
            
            
        }
        
    }
    
    void TimGraph::BuildHypergraphKPT(int64 R){
        
        // used for kpt estimation
        hyperG_temp.clear();
        for(int i = 0; i < n; i++)
            hyperG_temp.push_back(std::vector<int>());
        
        hyperGT_temp.clear();
        
        while((int)hyperGT_temp.size() <= R) // applies both to isTemp true and false
            hyperGT_temp.push_back(std::vector<int>());
        
        for(int i = 0; i < R; i++)
            BuildHypergraphNode(sfmt_genrand_uint32(&sfmtSeed)%n, i, true, hyperGT_temp);
        /*
         / / new cpp
         // 		for(int i = 0; i < R; i++){
         // 			for(int t:hyperGT_temp[i]) {
         // 				hyperG_temp[t].push_back(i);
         // 			}
         // 		}
         */
        
        // oldCPP
        for(int i = 0; i < R; i++){
            for(int j = 0; j < (int) hyperGT_temp[i].size(); j++) {
                int t = hyperGT_temp[i][j];
                hyperG_temp[t].push_back(i);
            }
        }
        
    }
    
    void TimGraph::estimateTheta() {
        
        // first estimateOPT via KptEstimation and refineKpt
        float ept, eps_prime;
        
        float lb=1/2.0;
        float c=0;
        int64 lastR=0;
        while(true){
            int loop= (6 * log(n)  +  6 * log(log(n)/ log(2)) )* 1/lb  ;
            c=0;
            lastR=loop;
            // 			IF_TRACE(int64 now=rdtsc());
            float sumMgTu=0;
            for(int i=0; i<loop; i++){
                int u=rand()%n;
                float MgTu=MgT(u);
                float pu=MgTu/m;
                sumMgTu+=MgTu;
                c+=1-pow((1-pu), kappa);
            }
            c/=loop;
            if(c>lb) break;
            lb /= 2;
        }
        BuildHypergraphKPT(lastR); // for refineKpt
        
        ept = (c * n) / 2; //KptEstimation
        // 		cout << "first step ept " << ept << endl;
        // estimateOPT - refineKpt kismi burasi
        // select kappa seeds from the hypergraph builded above : coming from the last iteration of KptEstimation
        BuildSeedSetTemp();
        eps_prime = 5.0 * std::pow(sqr(epsilon) / (kappa * 1.0), 1.0/3.0); //epsilon' used in RefineKPT
        int64 R = (8 + 2 * eps_prime) * ( n * log(n) +  n * log(2)  ) / ( eps_prime * eps_prime * ept)/4;
        BuildHypergraphKPT(R);
        
        std::set<int> s;
        
        for(set<int>::iterator t = seedSetTemp.begin(); t != seedSetTemp.end(); t++) {
            int u = *t;
            for(int j = 0; j < (int) hyperG_temp[u].size(); j++) {
                s.insert(hyperG_temp[u][j]);
            }
        }
        
        ept = (( float)(n * s.size())) / R;
        ept = ept / (1 + eps_prime);
        
        theta = (8 + 2 * epsilon) * ( n * log(n) + n * log(2) +  n * logcnk(n, kappa) ) / ( epsilon * epsilon * ept);
    }
    
    // will be used just to create the kappa seeds needed in RefineKPT
    void TimGraph::BuildSeedSetTemp() {
        std::vector<int> degree;
        std::vector<int> visit_local(hyperGT_temp.size());
        seedSetTemp.clear();
        
        for(int i = 0; i < n; i++) 	{
            degree.push_back(hyperG_temp[i].size());
        }
        
        for(int i = 0; i < kappa; i++){
            // new cpp
            // 			auto t = max_element(degree.begin(), degree.end());
            // 			int id = t - degree.begin();
            // oldCPP
            
            // ctr buraya gelebilir
            std::vector<int>::iterator t = max_element(degree.begin(), degree.end());
            int id = std::distance(degree.begin(), t);
            
            seedSetTemp.insert(id);
            degree[id]=0;
            // new cpp
            // 			for(int t : hyperG_temp[id]){
            // 				if(!visit_local[t]){
            // 					visit_local[t]=true;
            // 					for(int item : hyperGT_temp[t]){
            // 						degree[item]--;
            // 					}
            // 				}
            // 			}
            
            // oldCPP
            for(int j = 0; j < (int) hyperG_temp[id].size(); j++){
                int t = hyperG_temp[id][j];
                if(!visit_local[t]){
                    visit_local[t]=true;
                    for(int z = 0; z < (int) hyperGT_temp[t].size(); z++){
                        int item = hyperGT_temp[t][z];
                        degree[item]--;
                    }
                }
            }
            
            
            
        }
    }
    
    float TimGraph::logcnk(int n, int k){
        float ans = 0;
        
        for(int i = n - k + 1; i <= n; i++){
            ans += log(i);
        }
        for(int i = 1; i <= k; i++){
            ans -= log(i);
        }
        return ans;
    }
    
    int TimGraph::BuildHypergraphNode(int uStart, int hyperiiid, bool addHyperEdge, std::vector< std::vector<int> > &hyperGT){
        int n_visit_edge=1;
        
        if(addHyperEdge) {
            // 			ASSERT((int)hyperGT.size() > hyperiiid);
            hyperGT[hyperiiid].push_back(uStart);
        }
        
        int n_visit_mark=0;
        
        q.clear();
        q.push_back(uStart);
        // 		ASSERT(n_visit_mark < n);
        
        
        visit_mark[n_visit_mark++] = uStart;
        visit[uStart] = true;
        
        while(!q.empty()) {
            int expand = q.front();
            q.pop_front();
            
            int i = expand;
            for(int j = 0; j < (int)graphT[i].size(); j++){
                //int u=expand;
                int v = graphT[i][j]; //parent of u in the original graph G
                n_visit_edge++;
                float randFloat= float(sfmt_genrand_uint32(&sfmtSeed))/ float(RAND_MAX)/2;
                if(randFloat > adv->probT[i][j])
                    continue;
                if(visit[v])
                    continue;
                if(!visit[v]) {
                    // 					ASSERT(n_visit_mark < n);
                    visit_mark[n_visit_mark++]=v;
                    visit[v]=true;
                }
                q.push_back(v);
                
                if(addHyperEdge) {
                    // 					ASSERT((int)hyperGT.size() > hyperiiid);
                    hyperGT[hyperiiid].push_back(v);
                }
            }
        }
        
        for(int i = 0; i < n_visit_mark; i++)
            visit[visit_mark[i]]=false;
        
        return n_visit_edge; // returns number of edges considered, i.e., width of the RR set created from uStart
    }
    
    float TimGraph::MgT(int u){
        // 		ASSERT(u>=0);
        // 		ASSERT(u<n);
        return ( float)BuildHypergraphNode(u, 0, false, hyperGT_temp);
    }
    
    float TimGraph::sqr( float t) {
        return t*t;
    }
    
    
}

#endif

