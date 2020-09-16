//
// GeneExpressionAntithetic_h
//
//  Created by Ankit Gupta on 15.09.20
//Example file for the closed-loop network where the antithetic integral feedback controller controls the gene-expression network

#ifndef GeneExpressionAntithetic_h
#define GeneExpressionAntithetic_h
#endif /* GeneExpressionAntithetic_h */





#define NUMBER_OF_SPECIES 4
#define NUMBER_OF_REACTIONS 7
#define K_FB 1
#define OUTPUT_SPECIES 1

#define OUTPUT_FOLDER "./Results/Antithetic"

#define THETA 5
#define MU 100
#define ETA 1000


int StoichiometryMatrix[NUMBER_OF_SPECIES][NUMBER_OF_REACTIONS];
double propensity[NUMBER_OF_REACTIONS];
double ParameterValues[NUMBER_OF_REACTIONS];
int initialstate[NUMBER_OF_SPECIES] = {0,0,0,0}; // SPECIFY INITIAL STATE HERE
double PropensitySum;



// Objective function


 
double ProportionalFeedback(double x){
    double reference = MU/THETA;
    double y = K_FB*( 3*reference - x);
    if(y < 0){
        return 0;
    }
    else{
        return y;
    }
}

double HillFeedback(double x){
    double reference = MU/THETA;
    return 4*K_FB*reference*reference/(reference + x);
}


// propensity function
double CalculatePropensities(int * state){
        double x = ( (double) state[1] );
    propensity[0] = ParameterValues[0]*(state[2]) + ProportionalFeedback(x);
    propensity[1] = ParameterValues[1]*(state[0]);
    propensity[2] = ParameterValues[2]*(state[0]);
    propensity[3] = ParameterValues[3]*(state[1]);
    
    propensity[4] = MU;
    propensity[5] = THETA*(state[1]);
    propensity[6] = ETA*(state[2])*(state[3]);
    
    double sum = 0;
    
    for(int i = 0; i < NUMBER_OF_REACTIONS ; i++){
        sum = sum + propensity[i];
    }
    
    return sum;
}





void SetParameterValues(){
    
    // set the parameter values
    ParameterValues[0] =  4;
    ParameterValues[1] = 2;
    ParameterValues[2] = 5;
    ParameterValues[3] = 1;
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
    
    
    StoichiometryMatrix[2][4] = 1;
     StoichiometryMatrix[3][5] = 1;
    StoichiometryMatrix[2][6] = -1;
    StoichiometryMatrix[3][6] = -1;
    SetParameterValues();
}
