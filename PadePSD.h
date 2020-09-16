//
//  SpectrumEstimation.h
//
//  Created by Ankit Gupta on 15.09.20.
//

#ifndef PadePSD_h
#define PadePSD_h


using namespace std;


int NextreactionIndex;

#define NUM_SVALUES_VALIDATION 4
#define ApproximationOrder 2
#define MAX_DEPTH 10
#define BIG_PRIME_NUMBER  8022017 // This is a large prime number that determines the size of the hash table. It can be chosen according to the available memory resources. The bigger the better!


int HashFunction(int * state){
    int h = 216091; // a prime number
    int c = 132049; // another prime number
    for (int i=0 ; i < NUMBER_OF_SPECIES; ++i){
        h = (h*state[i] + c)%BIG_PRIME_NUMBER;
    }
    
    if(h < 0){
        h = h + BIG_PRIME_NUMBER;
    }
    
    return h;
}

// select the S - values for direct G function estimates. This can be used to compute the validation scores
double SValues[NUM_SVALUES_VALIDATION];
double SetSValues(double * SValues){
    double sum = 0;
    for(int i = 0; i < NUM_SVALUES_VALIDATION; i++){
        SValues[i] = 0.5*(i+1);
        sum = sum + SValues[i];
    }
    return sum;
}


double FindNextReactionSSA(int * state, double svalue_sum = 0){
    double TimeIncrement;
    double PropensitySum = CalculatePropensities(state) + svalue_sum;
    if(PropensitySum == 0){
        NextreactionIndex = -1;
        return DBL_MAX; // very big number
    }
    else{
        TimeIncrement = -log(Random(generator))/PropensitySum;
        double sum = 0;
        double unif = (Random(generator))*PropensitySum;
        
        for(int i = 0; i < (NUMBER_OF_REACTIONS + NUM_SVALUES_VALIDATION); i++){
            if(i < NUMBER_OF_REACTIONS){
                sum = sum + propensity[i];
            }
            else {
                sum = sum + SValues[(i-NUMBER_OF_REACTIONS )];
            }
            if(unif < sum){
                NextreactionIndex = i;
                break;
            }
        }
        return TimeIncrement;
    }
}

void GenerateSingleSSATrajectory(int * initialstate, int * finalstate, double initialtime, double finaltime){
    
    double TimeIncrement;
    double timer = initialtime;
    for(int i = 0; i < NUMBER_OF_SPECIES ; i++){
              finalstate[i] = initialstate[i];
      }
      while(1) {
          TimeIncrement = FindNextReactionSSA(finalstate);
          timer = timer + TimeIncrement;
          
          if(timer > finaltime){
              timer = finaltime;
              break;
          }
          else if (NextreactionIndex >=0){
              for(int i = 0 ; i < NUMBER_OF_SPECIES ; i++){
                  finalstate[i] = finalstate[i] + StoichiometryMatrix[i][NextreactionIndex];
              }
          }
      }
}


void GenerateSingleSSATrajectoryAugmented(int * initialstate, int * finalstate, int * finalAuxstate,double * Svalues,double svalue_sum, double initialtime, double finaltime){
    
    double TimeIncrement;
    double timer = initialtime;
    for(int i = 0; i < NUMBER_OF_SPECIES ; i++){
              finalstate[i] = initialstate[i];
    }
    for(int i = 0; i < NUM_SVALUES_VALIDATION ; i++){
              finalAuxstate[i] = finalstate[OUTPUT_SPECIES];
    }
      while(1) {
          TimeIncrement = FindNextReactionSSA(finalstate,svalue_sum);
          timer = timer + TimeIncrement;
          
          if(timer > finaltime){
              timer = finaltime;
              break;
          }
          else if (NextreactionIndex >=0 && NextreactionIndex < NUMBER_OF_REACTIONS){
              for(int i = 0 ; i < NUMBER_OF_SPECIES ; i++){
                  finalstate[i] = finalstate[i] + StoichiometryMatrix[i][NextreactionIndex];
              }
          }
          else if (NextreactionIndex >= NUMBER_OF_REACTIONS){
              finalAuxstate[(NextreactionIndex - NUMBER_OF_REACTIONS)] = finalstate[OUTPUT_SPECIES];
          }
      }
}

void GenerateSampledSSATrajectory(std::string filename, int * initialstate, double cutofftime, double finaltime, double Delta_t){

    std::ofstream SSAData;
    SSAData.open(filename);
    
    if (!SSAData.is_open()) {
        std::cout<<"File could not be created!"<<std::endl;
        exit(1);
    }
    int counter = 0;
    int state[NUMBER_OF_SPECIES];
    GenerateSingleSSATrajectory(initialstate,state,0,cutofftime);
    while( counter*Delta_t  <  (finaltime - cutofftime) ){
        for(int i = 0; i < NUMBER_OF_SPECIES ; i++){
            initialstate[i] = state[i];
        }
        GenerateSingleSSATrajectory(initialstate,state,counter*Delta_t + cutofftime,(counter+1)*Delta_t + cutofftime);
        counter++;
        SSAData<< (counter*Delta_t + cutofftime)<<" "<<state[OUTPUT_SPECIES]<<endl;
    }
   SSAData.close();
}
bool AreStatesEqual(int * a, int * b){
    for(int i = 0; i < NUMBER_OF_SPECIES; i++){
        if(a[i]!=b[i]){
            return false;
        }
    }
    return true;
}


struct NodeInfo{
    int storedstate[NUMBER_OF_SPECIES];
    bool GeneratorAvailability[2*ApproximationOrder];
    double GeneratorValues[2*ApproximationOrder];
    bool statefilled;
    int nextaddress;
    NodeInfo(){ //define the constructor
        statefilled = false;
    }
};
NodeInfo hashtable[2*BIG_PRIME_NUMBER];
int CurrentSpillOver = 0;

void AddState(int * state, int h){
        hashtable[h].statefilled = true;
        hashtable[h].nextaddress = -1;
        for(int j = 0; j < NUMBER_OF_SPECIES ; j++){
            hashtable[h].storedstate[j] = state[j];
        }
    
    for(int i=0; i < 2*ApproximationOrder; i++){
        hashtable[h].GeneratorAvailability[i] = false;
        hashtable[h].GeneratorValues[i] = -1;
    }
}
void PrintState(int * state){
        
    for(int j = 0; j < NUMBER_OF_SPECIES ; j++){
               cout<<state[j]<<" ";
           }
    
    cout<<endl;
}

int LookForStateInHashTable(int * state){ //looks for right location of the state, if not found it adds it
     int h0 = HashFunction(state);
    int address;
    int h=h0;
    if(!hashtable[h].statefilled){
         AddState(state,h);
    }
    else{
        int counter = 0;
         
        while(counter < MAX_DEPTH){
                if(AreStatesEqual(state,hashtable[h].storedstate)){
                    return h;
                }
                else if(hashtable[h].nextaddress>= 0){
                    h = hashtable[h].nextaddress + BIG_PRIME_NUMBER;
                }
                else{
                    break;
                }
            counter++;
            }
        
        if(hashtable[h].nextaddress >= 0 || (CurrentSpillOver >= BIG_PRIME_NUMBER) ){
            h = h0;
            address = hashtable[h].nextaddress;
            AddState(state,h);
            hashtable[h].nextaddress = address;
        }
        else {
            hashtable[h].nextaddress = CurrentSpillOver;
            h = BIG_PRIME_NUMBER + CurrentSpillOver;
            CurrentSpillOver++;
            AddState(state,h);
        }
    }
    return h;
}

double ComputeRecursiveGenerator(int * state, int n){ //A^n f(x) recursive computation
    
    if(n == 0){
        return (double) state[OUTPUT_SPECIES];
    }
  
    int h = LookForStateInHashTable(state);
  
    double value = 0;
    
    if(hashtable[h].GeneratorAvailability[n]){
         value = hashtable[h].GeneratorValues[n];
    }
    else{
        double newpropensity[NUMBER_OF_REACTIONS];
        double sum;

        sum = CalculatePropensities(state);
        for(int i = 0; i < NUMBER_OF_REACTIONS ; i++){
            newpropensity[i] = propensity[i];
        }
        value = - sum*ComputeRecursiveGenerator(state,n-1);

        for(int i = 0; i < NUMBER_OF_REACTIONS ; i++){
            if (newpropensity[i] > 0){
                int newstate[NUMBER_OF_SPECIES]= {0};
                for(int j = 0; j < NUMBER_OF_SPECIES ; j++){
                    newstate[j] = state[j] + StoichiometryMatrix[j][i];
                }
                    value = value + newpropensity[i]*ComputeRecursiveGenerator(newstate,n-1);
                }
        }
        hashtable[h].GeneratorAvailability[n] = true;
        hashtable[h].GeneratorValues[n] = value;
        
    }
    return value;
}
double ComputeGammaFunction(int * state, int n1, int n2){
    double newpropensity[NUMBER_OF_REACTIONS];
    double sum = CalculatePropensities(state);
    
    for(int i = 0; i < NUMBER_OF_REACTIONS ; i++){
        newpropensity[i] = propensity[i];
    }
    double value = 0;
    double fnbase1  = ComputeRecursiveGenerator(state,n1);
    double fnbase2 = fnbase1;
    if(n1!=n2){
        fnbase2 =  ComputeRecursiveGenerator(state,n2);
    }
    int newstate[NUMBER_OF_SPECIES]= {0};
    
    for(int i = 0; i < NUMBER_OF_REACTIONS ; i++){
        if (newpropensity[i] > 0){
            for(int j = 0; j < NUMBER_OF_SPECIES ; j++){
                newstate[j] = state[j] + StoichiometryMatrix[j][i];
            }
            if(n1!=n2){
                value = value + newpropensity[i]*(ComputeRecursiveGenerator(newstate,n1) - fnbase1)*(ComputeRecursiveGenerator(newstate,n2) - fnbase2);
            }
            else{
                value = value + newpropensity[i]*pow( ComputeRecursiveGenerator(newstate,n1) - fnbase1,2);
            }
        }
    }
    return value;
}

double nchoosek(int a , int b){

    if(a < b){
        return 0;
    }
    else if(a==b || b == 0 ) {
        return 1;
    }
    else{
        return nchoosek(a-1,b-1) + nchoosek(a-1,b);
    }
}

void  SingleTrajectory_PadeDerivatives(std::string trajectory_filename, std::string DRM_filename, std::string IGM_filename, int * initialstate, double cutofftime, double finaltime, double Delta_t){

    std::ofstream IGM;
    std::ofstream DRM;
    std::ofstream SSAData;
    
    SSAData.open(trajectory_filename);
    IGM.open(IGM_filename);
    DRM.open(DRM_filename);
    if (!IGM.is_open() || !DRM.is_open() || !SSAData.is_open() ) {
        std::cout<<"File could not be created!"<<std::endl;
        exit(1);
    }
    
    int counter = 1;
    int TrajectoryCounter = 0;
    double TimeIntervalPrint = 10;
    double TimeIncrement;
    double timer = cutofftime;
    int state[NUMBER_OF_SPECIES];
    int Auxstate[NUM_SVALUES_VALIDATION];
    
    bool stop_flag = false;
    double incr;
    double sum_mean = 0;
    double sum_fvalue[2*ApproximationOrder] = {0};
    double sum_resolvents[NUM_SVALUES_VALIDATION] = {0};
    double fvalue[(2*ApproximationOrder-1)];
    int m ;
 
    double svalue_sum = SetSValues(SValues);
    
    GenerateSingleSSATrajectoryAugmented(initialstate,state,Auxstate,SValues,svalue_sum,0,cutofftime);
    cout<<"Cutoff period ended!"<<endl;
    //SSAData<<"Time"<<"   Output"<<endl;
    while(!stop_flag) {
     
        TimeIncrement = FindNextReactionSSA(state,svalue_sum);
        if(timer > (finaltime - TimeIncrement)){
            TimeIncrement = finaltime - timer;
            stop_flag = true;
        }
        
        while( (timer > cutofftime) && ((TrajectoryCounter*Delta_t + cutofftime) < (timer + TimeIncrement) ) ){
            SSAData << (TrajectoryCounter*Delta_t + cutofftime)<<" "<<state[OUTPUT_SPECIES];
            for( int n = 0; n < 2*ApproximationOrder; n++){
                SSAData<<" "<<sum_fvalue[n]/(timer - cutofftime);
            }
            SSAData<<endl;
            TrajectoryCounter++;
        }
     
        for( int n = 0; n < (2*ApproximationOrder-1); n++){
            fvalue[n] = ComputeRecursiveGenerator(state,n);
         }
     
        sum_fvalue[0] = sum_fvalue[0] + TimeIncrement*pow(fvalue[0],2);
        for( int n = 1; n < 2*ApproximationOrder; n++){
            incr = 0;
            if( (n % 2) == 0){
                m = (int) n/2;
                for(int k = 0; k < m; k++){
                    if(k > 0){
                        incr = incr + nchoosek(n,k)*fvalue[k]*fvalue[n-k];
                    }
                    incr = incr + nchoosek(n-1,k)*ComputeGammaFunction(state,k,n-1-k);
                }
                incr = incr + 0.5*nchoosek(n,m)*pow(fvalue[m],2);
            }
            else {
                m = (int) (n-1)/2;
                for(int k = 1; k <=m; k++){
                    incr = incr + nchoosek(n,k)*fvalue[k]*fvalue[n-k];
                    incr = incr + nchoosek(n-1,k-1)*ComputeGammaFunction(state,k-1,n-k);
                }
                incr = incr + 0.5*nchoosek(n-1,m)*ComputeGammaFunction(state,m,m);
            }
            sum_fvalue[n] = sum_fvalue[n] - incr*TimeIncrement;
        }
        
        for( int n = 0; n < NUM_SVALUES_VALIDATION; n++){
            sum_resolvents[n] = sum_resolvents[n] + Auxstate[n]*TimeIncrement*fvalue[0];
        }
        
        sum_mean = sum_mean + TimeIncrement*(fvalue[0]);
        timer = timer + TimeIncrement;
        
        if(timer > counter*TimeIntervalPrint + cutofftime){
            cout<<" Reached Time "<<counter*TimeIntervalPrint<<endl;
            counter++;
        }
        if (NextreactionIndex >=0 && NextreactionIndex < NUMBER_OF_REACTIONS){
            for(int i = 0 ; i < NUMBER_OF_SPECIES ; i++){
                state[i] = state[i] + StoichiometryMatrix[i][NextreactionIndex];
            }
        }
        else if (NextreactionIndex >= NUMBER_OF_REACTIONS){
            Auxstate[(NextreactionIndex - NUMBER_OF_REACTIONS)] = state[OUTPUT_SPECIES];
        }
    }
    double TimeIntervalLength = (finaltime - cutofftime);
    sum_mean =  sum_mean/TimeIntervalLength;
    cout<<"Mean Value "<<sum_mean<<endl;

    IGM<<"Order Estimate"<<std::endl;
    for( int n = 0; n < 2*ApproximationOrder; n++){
        sum_fvalue[n] =  sum_fvalue[n]/TimeIntervalLength;
        if(n==0){
            sum_fvalue[n] = sum_fvalue[n] - pow(sum_mean,2);
        }
        IGM<<n<<" "<<sum_fvalue[n]<<endl;
    }
    //cout<<"Variance Value "<<sum_fvalue[0]<<endl;
    //cout<<"Fano factor "<<sum_fvalue[0]/sum_mean<<endl;
   IGM.close();
    
    DRM<<"Svalues Resolvent_Estimate"<<std::endl;
      for( int n = 0; n < NUM_SVALUES_VALIDATION; n++){
          sum_resolvents[n] =  sum_resolvents[n]/TimeIntervalLength;
          sum_resolvents[n] = (sum_resolvents[n] - pow(sum_mean,2))/SValues[n];
          DRM<<SValues[n]<<" "<<sum_resolvents[n]<<endl;
      }
    DRM.close();
    SSAData.close();
}




#endif /* SpectrumEstimation_h */
