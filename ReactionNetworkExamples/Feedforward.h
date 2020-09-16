
//  Feedforward.h
//  Created by Ankit Gupta on 15.09.20.
//Example file for the Feedforward Network


#ifndef Feedforward_h
#define Feedforward_h


#endif /* Feedforward_h */

#define gamma_c 0.5
#define gamma_o 0.3
#define k_c 0.25
#define k_o 2
#define beta_0 50
#define beta_ff 4
#define I_0 3

#define NUMBER_OF_SPECIES 2
#define NUMBER_OF_REACTIONS 4
#define OUTPUT_SPECIES 1
#define OUTPUT_FOLDER "./Results/Feedforward"

int StoichiometryMatrix[NUMBER_OF_SPECIES][NUMBER_OF_REACTIONS];
double propensity[NUMBER_OF_REACTIONS];
int initialstate[NUMBER_OF_SPECIES] = {0}; // SPECIFY INITIAL STATE HERE




// propensity function
double CalculatePropensities(int * state){
    propensity[0] = k_c*I_0;
    propensity[1] = k_o*I_0 +  beta_0 - beta_ff*(state[0]);
    if(propensity[1] < 0){
         propensity[1] = 0;
     }
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

