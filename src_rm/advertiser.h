#ifndef _ADVERTISER_H
#define _ADVERTISER_H

namespace _Cide{

    class advertiser{
        
    public:
        
        std::vector< std::vector< float> > probT;
        int advertiserID;
        float *gamma;
        std::vector<float> seedUserCosts; // advertiser-specific
        float maxCost; // to be used in the timgraph object
        float minCost; // to be used in the timgraph object 
        
        float currentPayment;
        float currentRevenue;
        float currentSeedCosts;
        
        float budget, cpe;
        
        advertiser(int id, int nrTopics) {
            this->advertiserID = id;
            this->gamma = new  float[nrTopics];
            
        }
        
        ~advertiser() {}
        
        void setItemDist( float *temp, int size) {
            for(int i = 0; i < size; i++)
                this->gamma[i] = temp[i];
        }
    
    
    };


}


#endif
