
//  Feedback.h
//  Created by Ankit Gupta on 15.09.20.
//Example file for the Negative Feedback Network

#ifndef Feedback_h
#define Feedback_h


#endif /*Feedback_h */


//define network specific parameters
#define gamma_c 1
#define gamma_o 0.5
#define k_o 2
#define beta_0  200
#define beta_fb 0.5
#define I_0 1

#define NUMBER_OF_SPECIES 2
#define NUMBER_OF_REACTIONS 4
#define OUTPUT_SPECIES 1
#define INPUT_SPECIES 0

#define OUTPUT_FOLDER "./Results/Feedback"


int StoichiometryMatrix[NUMBER_OF_SPECIES][NUMBER_OF_REACTIONS];
double propensity[NUMBER_OF_REACTIONS];
int initialstate[NUMBER_OF_SPECIES] = {0}; // SPECIFY INITIAL STATE HERE






// propensity function
double CalculatePropensities(int * state){
    propensity[0] = I_0*(beta_0 - beta_fb*(state[1]));
    if(propensity[0] < 0){
        propensity[0] = 0;
    }
    
    
    propensity[1] = k_o*(state[0]);
    propensity[2] = gamma_c*(state[0]);
    propensity[3] = gamma_o*(state[1]);
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
    
    // set the stoichiometry matrix. Only specify the non-zero entries
    StoichiometryMatrix[0][0] = 1;
    StoichiometryMatrix[0][2] = -1;
    StoichiometryMatrix[1][1] = 1;
    StoichiometryMatrix[1][3] = -1;
    
}



