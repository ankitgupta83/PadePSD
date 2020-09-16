
//  SelfRegulatory.h
//  Created by Ankit Gupta on 15.09.20.
//Example file for the Self Regulatory gene expression model consider in Figure 1

#ifndef SelfRegulatory_h
#define SelfRegulatory_h
#endif /* SelfRegulatory_h */





#define NUMBER_OF_SPECIES 1
#define NUMBER_OF_REACTIONS 2
#define OUTPUT_SPECIES 0
#define HILL_FACTOR 1 // 0 for no regulation, 1 for regulation

#define OUTPUT_FOLDER "./Results/SelfRegulatory/WithRegulation"





int StoichiometryMatrix[NUMBER_OF_SPECIES][NUMBER_OF_REACTIONS];
double propensity[NUMBER_OF_REACTIONS];
double ParameterValues[NUMBER_OF_REACTIONS];
int initialstate[NUMBER_OF_SPECIES] = {0}; // SPECIFY INITIAL STATE HERE
double PropensitySum;




// propensity functions
double CalculatePropensities(int * state){
        double x = ( (double) state[1] );
    propensity[0] = (2-HILL_FACTOR)*pow(ParameterValues[0],HILL_FACTOR)/(pow(ParameterValues[0],HILL_FACTOR) + pow(state[0],HILL_FACTOR)); ;
    propensity[1] = ParameterValues[1]*(state[0]);
    double sum = 0;
    
    for(int i = 0; i < NUMBER_OF_REACTIONS ; i++){
        sum = sum + propensity[i];
    }
    
    return sum;
}

void SetParameterValues(){
    
    // set the parameter values
    ParameterValues[0] =  10;
    ParameterValues[1] = 2;
}
void SetStoichiometryMatrix(){
    for(int i = 0 ; i < NUMBER_OF_SPECIES ; i++){
        for(int j = 0 ; j < NUMBER_OF_REACTIONS ; j++ ){
            StoichiometryMatrix[i][j] = 0;
        }
    }
    
    // set the stoichiometry matrix. Only specify the non-zero entries
    StoichiometryMatrix[0][0] = 1;
    StoichiometryMatrix[0][1] = -1;
    SetParameterValues();
}
