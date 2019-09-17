#ifndef _TIMGRAPH_H_
#define _TIMGRAPH_H_

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#include <vector>
#include <deque>
#include <utility>
#include "sfmt/SFMT.h"
#include "advertiser.h"
#include "utils.h"
#include <unordered_set>
#include <map>

#define IF_TRACE(args) ;

namespace _Cide {
	
	
	class TimGraph {
		
	public: 
		
		//std::vector< std::vector< float> > probT; // keep the item specific influence probs in this pbject, in transposed graph style
		//std::vector<float> nodeCTRs;
		
		// maintain RR sets tru execution
		std::vector< std::vector<int> > hyperG_adv;
		std::vector< std::vector<int> > hyperGT_adv;
		
		// keep the temp versions to be used in the kpt estimations
		std::vector< std::vector<int> > hyperG_temp;
		std::vector< std::vector<int> > hyperGT_temp;		
		std::set<int> seedSetTemp; // temporary seed set created for Kpt estimation purposes	
		
		// acts local
		std::vector<bool> visit; 
		std::vector<int> visit_mark; 
		
		std::deque<int> q;
		sfmt_t sfmtSeed;	
		
		// needs to be actively updated
		std::vector<int> seedSet; // keep like this to be able to update in this order
		std::vector<int> nrRRCovered; //keep this with number of covered (not m.g.) to be able to update later
		
		std::vector<int> hyper_degree; // 0 to n, her node icin var - bu RR set yaratilirken doldurulmali, hyperG_adv'den dolduruluyor (hyperG_adv[i].size())
		std::vector<bool> isCovered; // her RR_id icin - bu da RR set yapilirken de
		
		advertiser *adv;
		int candidateNode;
		float candidateMgRev; // candidate mg rev
		int candidateNrRR; // candidate nr of covered 
		int64 theta;
		int64 theta_old;
		int kappa;		
		float epsilon;
        float lambda;
		int n, m;
        int windowSize;
        
        // optim denemesi
        multimap<float,int> criterQueue; // the float keys differ wrt cost-sensitive or not
        set<int> usersExamined;
        float mmg_temp;
		
		TimGraph(advertiser *adv,  float eps, int nrNodes, int nrEdges, string greedyKriter, int windowSize);
		~TimGraph(void);		
		
		void doInitialGeneration();
		void GenerateRRSets();
		void findBestCANode(infPair &bestCandidate);
        void findBestCSNode(infPair &bestCandidate);
		void assignBestNode();
		void updateEstimates();
		
		int BuildHypergraphNode(int uStart, int hyperiiid, bool addHyperEdge, std::vector< std::vector<int> > &hyperGT);
        
        // yeni
        std::vector<float> seedMgRevs;
        float max_cost; // needed for the latent seed set size estimation procedure (both for CA and CS)
        
        string criterGreedy; 
		
		void estimateTheta();
		void BuildHypergraphKPT(int64 R);
		void BuildSeedSetTemp();
		float MgT(int u);
		float logcnk(int n, int k);
		float sqr( float t);
        int vec_intersect(const std::vector<int>& A_, const std::vector<int>& B_);
		
		uint64 rdtsc(void)
		{
			unsigned a, d;
			//asm("cpuid");
			asm volatile("rdtsc" : "=a" (a), "=d" (d));
			return (((uint64)a) | (((uint64)d) << 32));
		}
		
	};
	
}


#endif
