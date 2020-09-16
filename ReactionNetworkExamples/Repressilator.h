
//  Repressilator.h
//  Created by Ankit Gupta on 15.09.20.
//Example file for the Repressilator model

#ifndef Repressilator_h
#define Repressilator_h


#endif /* Repressilator_h */
#define lambda 60
#define K1 5
#define K2 10
#define K3 10
#define Mean_No 10
#define Mean_Nt 0
#define Mean_B_Inv 1 // enter reciprocal of the mean of the burst size

#define NUMBER_OF_SPECIES 5  // three repressor proteins, N_o and N_t
#define NUMBER_OF_REACTIONS 10
#define OUTPUT_SPECIES 1
#define HILL_FACTOR 2 // This models the cooperativity
#define OUTPUT_FOLDER "./Results/Repressilator"



std::geometric_distribution<int> GeomRand(Mean_B_Inv);
int set_geometric_rnd(){
  return (1 + GeomRand(generator));
}

 

int StoichiometryMatrix[NUMBER_OF_SPECIES][NUMBER_OF_REACTIONS];
double propensity[NUMBER_OF_REACTIONS];
int initialstate[NUMBER_OF_SPECIES] = {0}; // SPECIFY INITIAL STATE HERE

int set_geometric_rnd();



double ComputeFreeProteins(int N_pro,double p_tot,double K,double hill_coeff){ // we solve the nonlinear equation with the bisection method
    //using bisection method
    double tol =0.01;
    double err;
    double p_free = 0;
    
    double left = 0;
    double right = p_tot;

    while(1){
      // std::cout<<p_free<<" "<<p_tot<<" "<<left<<" "<<right<<" "<<N_pro<<std::endl;     
        p_free = (left+right)/2;
        err = p_tot - p_free - 2*N_pro*pow(p_free,hill_coeff)/( pow(p_free,hill_coeff) + pow(K,hill_coeff) );
        if (err < -tol){
            right = p_free;
        }
        else if (err > tol){
            left = p_free;
        }
        else{
            return p_free;
        }
    }
}



// propensity function
double CalculatePropensities(int * state){
    
    
    StoichiometryMatrix[0][0] = set_geometric_rnd();
    StoichiometryMatrix[1][1] = set_geometric_rnd();
    StoichiometryMatrix[2][2] = set_geometric_rnd();
    
    double p_tot_1 = (double) state[0];
    double p_tot_2 = (double) state[1];
    double p_tot_3 = (double) state[2];
    
    // compute p_free_i
    double p_free_1 = ComputeFreeProteins(state[3]+state[4],p_tot_1,K1,HILL_FACTOR);
    double p_free_2 = ComputeFreeProteins(state[3],p_tot_2,K2,HILL_FACTOR);
    double p_free_3 = ComputeFreeProteins(state[3],p_tot_3,K3,HILL_FACTOR);
    
    //regulated production of repressor proteins
    propensity[0] = lambda*Mean_No*(pow(K1,HILL_FACTOR)/( pow(K1,HILL_FACTOR) + pow(p_free_3,HILL_FACTOR)));
    propensity[1] = lambda*Mean_No*(pow(K2,HILL_FACTOR)/( pow(K2,HILL_FACTOR) + pow(p_free_1,HILL_FACTOR)));
    propensity[2] = lambda*Mean_No*(pow(K3,HILL_FACTOR)/( pow(K3,HILL_FACTOR) + pow(p_free_2,HILL_FACTOR)));
    
    //degradation/dilution of repressor proteins
    propensity[3] = 1*p_tot_1;
    propensity[4] = 1*p_tot_2;
    propensity[5] = 1*p_tot_3;
    
    
    //plasmid production
    propensity[6] = Mean_No;
    propensity[7] = Mean_Nt;
    
    //plasmid degradation
    propensity[8] = state[3];
    propensity[9] = state[4];

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
    StoichiometryMatrix[1][1] = 1;
    StoichiometryMatrix[2][2] = 1;
    
    StoichiometryMatrix[0][3] = -1;
    StoichiometryMatrix[1][4] = -1;
    StoichiometryMatrix[2][5] = -1;
    
    StoichiometryMatrix[3][6] = 1;
    StoichiometryMatrix[4][7] = 1;
    
    StoichiometryMatrix[3][8] = -1;
    StoichiometryMatrix[4][9] = -1;
    
}
