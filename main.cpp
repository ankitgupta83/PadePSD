//
//  main.cpp
//  Padé PSD method
//  Created by Ankit Gupta on 15.09.20.

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <random>
#include <vector>
#include <math.h>
#include <string>
#include <cfloat>


#define NUM_TRAJ_FOR_NUMERICAL_PSD_AVERAGING 10 //this is for generating trajectories for the averaged periodogram method


//select uniform random number generator
std::random_device rd;
std::mt19937 generator(rd());
std::uniform_real_distribution<double> Random(0.0,1.0);


#include "./ReactionNetworkExamples/Feedforward.h"
//#include "./ReactionNetworkExamples/Feedback.h"
//#include "./ReactionNetworkExamples/Repressilator.h"
//#include "./ReactionNetworkExamples/Splicing.h"
//#include "./ReactionNetworkExamples/GeneExpressionAntithetic.h"
//#include "./ReactionNetworkExamples/SelfRegulatory.h"
#include "./PadePSD.h"
using namespace std;




int main(int argc, const char * argv[]) {
    SetStoichiometryMatrix();
    // set simulation parameters
    double TotalTime = 10000;
    double CutoffTime = 100;
    double Delta_t = 0.1;
    double FinalTime = CutoffTime + TotalTime;
    int recordlength = (int) ( TotalTime/Delta_t );
    int K_up = (int)(recordlength/2);
    int NumTraj = NUM_TRAJ_FOR_NUMERICAL_PSD_AVERAGING;
    
    

    string Pade_filename, DirectG_Estimates_filename, trajectory_filename, autocov_filename;
    // declare filenames
    Pade_filename = OUTPUT_FOLDER;
    DirectG_Estimates_filename = OUTPUT_FOLDER;
    trajectory_filename = OUTPUT_FOLDER;
    autocov_filename = OUTPUT_FOLDER;
    Pade_filename =  Pade_filename + "/PadeDerivatives.txt";
    DirectG_Estimates_filename =  DirectG_Estimates_filename + "/DirectG_Estimate.txt";
    trajectory_filename = trajectory_filename +  "/StochasticTrajectory.txt";
    autocov_filename = autocov_filename + "/AutoCovariance.txt";

    
    

    unsigned long start_s = clock();
    SingleTrajectory_PadeDerivatives(trajectory_filename, DirectG_Estimates_filename, Pade_filename, initialstate, CutoffTime,FinalTime, Delta_t);
    
    // uncomment the next for-loop for generating sampled- SSA trajectories for the Averaged Periodogram PSD estimation method
    
    /*
    for(int i = 0; i < NumTraj; i++){
        cout<<"Generating trajectory "<<(i+1)<<endl;
        string newtrajectoryfilename = OUTPUT_FOLDER;
        newtrajectoryfilename = newtrajectoryfilename + "/StochasticTrajectory" + to_string(i+1) + ".txt";
        GenerateSampledSSATrajectory(newtrajectoryfilename,initialstate,CutoffTime,FinalTime,Delta_t);
    }
    */
    return 0;
}
