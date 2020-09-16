// GeneExpressionAntithetic_h
//
//  Created by Ankit Gupta on 15.09.20.
//Example file for the RNA Splicing network
#ifndef Splicing_h
#define Splicing_h


#endif /* Splicing_h */



#define gamma 0.5
#define beta 2
#define k_on 1
#define k_off 3
#define alpha_on 3
#define alpha_off 0.2

#define NUMBER_OF_SPECIES 4
#define NUMBER_OF_REACTIONS 5
#define OUTPUT_SPECIES 3
#define OUTPUT_FOLDER "./Results/Splicing"

int StoichiometryMatrix[NUMBER_OF_SPECIES][NUMBER_OF_REACTIONS];
double propensity[NUMBER_OF_REACTIONS];
int initialstate[NUMBER_OF_SPECIES] = {1,0,0,0}; // SPECIFY INITIAL STATE HERE



// propensity function
double CalculatePropensities(int * state){
    
    

    
    propensity[0] = k_on*( (double) state[0]); // gene switching on
    propensity[1] = k_off*( (double) state[1]); // gene switching off
    propensity[2] = alpha_off + (alpha_on  - alpha_off)*((double) state[1]); // transcription
    propensity[3] = beta*((double) state[2]); //splicing
    propensity[4] = gamma*((double) state[3]); // degradation
    
    double sum = 0;
    for(int i = 0; i < NUMBER_OF_REACTIONS ; i++){
        sum = sum + propensity[i];
    }
    return sum;
}

void SetStoichiometryMatrix(){
    for(int i = 0 ; i < NUMBER_OF_SPECIES ; i++){
        for(int j = 0 ; j < NUMBER_OF_REACTIONS ; j++ ){
            StoichiometryMatrix[i][j] = 0;
        }
    }

    StoichiometryMatrix[1][0] = 1;
    StoichiometryMatrix[0][0] = -1; // gene switching on
    
    StoichiometryMatrix[1][1] = -1;
    StoichiometryMatrix[0][1] = 1; // gene switching off
     
    StoichiometryMatrix[2][2] = 1; // transcription
    
    StoichiometryMatrix[2][3] = -1;
    StoichiometryMatrix[3][3] = 1;
    
    StoichiometryMatrix[3][4] = -1;
    
}
