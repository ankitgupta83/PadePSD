# PadePSD
Implementation for the Pade PSD method developed in the paper 'Frequency Spectra and the Color of Cellular Noise' by Ankit Gupta and Mustafa Khammash

This repository contains the following files:

1. PadePSD_Plot.m -  the Matlab script which uses estimates of the Padé derivatives to  compute the PSD. It also uses direct estimates of the function G(s) to compute the validation score for the estimate. Validation scores close to 1 suggest that the estimated PSD is accurate.

2. main.cpp - this is the main C++ file, where one can select the reaction network (specified via a header file), choose simulation parameters (cutoff time, final time etc.) and can also generate discrete-sampled trajectory (with time discretisation specified by Delta_t) for estimating PSD with the averaged periodogram approach.

3. PadePSD.h - this is the header file which contains the implementation of Pade PSD method as described in the Supplement of the paper. One can choose the approximation order and the Hash table parameters at the top. 

4. ReactionNetworkExamples -  this is a folder with header files describing all the examples in the paper. Each file specifies the propensity functions and the stoichiometry vectors for all the reactions. It also specifies the output folder where the text files containing the estimates of the Padé derivatives and direct G-function estimates are kept. These text files are read by the Matlab script PadePSD_Plot.m to produce the PSD plot and compute the validation score. Make sure that the specified output folder exists!

If you have any questions regarding the code, please contact Ankit Gupta at ankit.gupta@bsse.ethz.ch.
