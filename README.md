# Simulation of evolution without mutation in coding sequences

This package aims at simulating evolution in absence of mutation in gene or promoter sequences. It is based on the work by [El Houdaigui et al., NAR 2019](https://doi.org/10.1093/nar/gkz300). 

## Getting Started

Clone the repository and run ourCode.py

## Inputs, Outputs & Parameters files :

All parameter files are in TwisTranscripT folder. 

When you run the code, the following parameters will be asked:

`path and name of the output csv file` :
Type in the output csv file (ex: out.csv) or press enter if you do not want to save the results

`Unit of length in base pairs that is deleted or inserted` :
Type in an integer number. Please use a multiple of the simulation discretization length (60) and do not exceed 300

`Probability for an evolutive event to be an inversion` :
Type in a float between 0 and 1

`Probability for an indel event to be an insertion`:
Type in a float between 0 and 1

`Number of generations` :
Type in an integer. 500 generations run about 45 minutes.

`Value of q` :
Type in a value for q. Under 1e-5 will refuse all negative mutations, and above 2e-3 will accept almost all of them.
