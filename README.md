# Chi-bar, Atomic and Molecular Electronegativity
Generation of the function nchi for molecular systems from Gaussian calculations

 X_space Version 2.0 October 2024

 Written by Dr. Stefano Racioppi, SUNY Buffalo. email: sraciopp@buffalo.edu

 This program generates a cube file containing the spatial distribution of X(r).
 following a quantum mechanical calculation of a molecule performed with the software Gaussian.
 Both closed and open-shell wavefunctions are OK.
 Requirements: numpy, cclib, argparse and ase >Python 3.0. 

 Please cite: https://doi.org/10.1002/chem.202103477

USAGE:
 use "python X_Space.py -i inputfile -o outputfile.cube" 

Script description:
 Example usage: python X_space.py -i file.log -f filelist.dat -n 6 -o outputFile.read_cube_data
Description  -f name of file with list of cube files that should be summed
             -i name of output file from the molecular calculation
             -n grid width
             -o name of outputfile
             
Further instructions are reported in the Supporting Information of https://doi.org/10.1002/chem.202103477
