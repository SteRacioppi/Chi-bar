import os
import numpy as np
from ase.io.cube import read_cube_data
import argparse
import cclib

# Argument parser setup
parser = argparse.ArgumentParser(description="Rescale cube files based on molecular orbital energies and sum them.")
parser.add_argument('-i', '--inputfile', default="file.log", help="Input file from molecular calculation")
parser.add_argument('-n', '--ngridpoints', type=int, default=6, help="Number of grid points")
parser.add_argument('-f', '--filelist', default="orbitalslist.dat", help="List of cube files to sum")
parser.add_argument('-o', '--output', default="output.cube", help="Output cube file")
args = parser.parse_args()

# File and parameter setup
filename = args.inputfile
fileList = args.filelist
gridpoints = args.ngridpoints
outfileName = args.output

# Parse molecular orbital data using cclib
cclibparser = cclib.io.ccopen(filename)
parsed_data = cclibparser.parse()

# Check if 'homos' is a list or an integer
if isinstance(parsed_data.homos, list):
    openshell = (len(parsed_data.homos) == 2)
else:
    openshell = False  # Default to closed-shell if it's not a list

# Extract and write molecular orbital energies to energies.dat
with open("energies.dat", "w") as output_data:
    for i in range(parsed_data.homos[0] + 1):
        print(parsed_data.moenergies[0][i] * -1, file=output_data)
    if openshell:
        for i in range(parsed_data.homos[1] + 1):
            print(parsed_data.moenergies[1][i] * -1, file=output_data)

# Load energies and file list
constantList = np.loadtxt("energies.dat")
filenameList = np.loadtxt(fileList, dtype='str')

# Main processing loop
data_sum = None
for counter, filename in enumerate(filenameList):
    multConstant = constantList[counter]  # Define multConstant at the beginning of the loop
    
    # Initialize a flag to track if corrections were made
    corrections_made = False
    
    while True:  # Keep checking until we either succeed or exhaust correction options
        try:
            # Attempt to read the cube data
            data, atoms = read_cube_data(filename)

            # If successful, break the loop
            break

        except ValueError as e:
            if "negative dimensions are not allowed" in str(e):
                # If the error is due to negative dimensions, adjust the first value in the third line of the header
                with open(filename, 'r') as f:
                    lines = f.readlines()

                # Change the first value of the third line to its absolute value
                atom_info = lines[2].split()
                atom_count = int(atom_info[0])
                if atom_count < 0:
                    atom_info[0] = str(abs(atom_count))  # Change to absolute value
                    lines[2] = " ".join(atom_info) + "\n"

                    # Rewrite the file with the corrected atom count
                    with open(filename, 'w') as f:
                        f.writelines(lines)

                    print(f"Negative atom count encountered in {filename}. I edited the cub file to fix this error.")
                    corrections_made = True  # Mark that a correction was made

            elif "cannot reshape array of size" in str(e):
                print(f"ValueError occurred while processing {filename}: {e}")
                print("It is probable that you have one more line in your header. I will try to fix it.")

                # Read the cube file to adjust the header
                with open(filename, 'r') as f:
                    lines = f.readlines()

                # Get the absolute atom count from the third line
                atom_info = lines[2].split()
                atom_count = abs(int(atom_info[0]))  # Atom count from third line

                # Calculate the expected number of header lines
                expected_header_lines = atom_count + 6
                
                # If there are more lines than expected, remove the next line
                if len(lines) > expected_header_lines:
                    print(f"Removing the 11th line from {filename} to fix the header. I edited the cub file to fix this error.")
                    del lines[expected_header_lines]  # Remove the next line after header

                    # Rewrite the adjusted cube file
                    with open(filename, 'w') as f:
                        f.writelines(lines)
                    
                    corrections_made = True  # Mark that a correction was made

            # If no corrections were made, break from the loop
            if not corrections_made:
                print(f"ValueError occurred while processing {filename}: {e}")
                break  # Exit the loop on unhandled errors

    # After all corrections, check if the file is valid
    if corrections_made:
        # Retry reading the cube data after correction
        try:
            data, atoms = read_cube_data(filename)
        except ValueError as retry_error:
            print(f"ValueError occurred while processing {filename} after correction: {retry_error}")
            continue  # Skip to the next file

    # If we successfully read the data (either directly or after corrections)
    if data is not None:
        # Calculate the squared modulus of each element in data and apply the multiplication
        sqMod = (data ** 2) * (2 if not openshell else 1)
        multsqMod = sqMod * multConstant  # Perform multiplication once here
        print(f"{multConstant:.3f} * {filename}")

        # Sum the rescaled cube data
        if counter == 0:
            data_sum = multsqMod
        else:
            data_sum += multsqMod

# Write the summed data to the output cube file
if data_sum is not None:
    with open(outfileName, "w") as fileOut:
        # Write header from the first cube file
        with open(filenameList[0], "r") as firstCube:
            for line in firstCube:
                if "E-" not in line:
                    fileOut.write(line)
                else:
                    break

        # Reshape and write the grid data
        nElements = data.shape[1] * data.shape[2]
        nFrames = data.shape[0]
        data_sum_reshape = data_sum.reshape(1, nFrames, nElements)

        for i in range(nFrames):
            frame = data_sum_reshape[0, i]
            for j in range(nElements):
                element = f'{frame[j]:1.9E}'
                fileOut.write(f"{element}\n" if (j + 1) % gridpoints == 0 or j == nElements - 1 else f"{element} ")

    # Modify the title in the output file
    with open(outfileName, "r") as fileOut:
        lines = fileOut.readlines()

    lines[1] = "NCHI\n"

    with open(outfileName, "w") as fileOut:
        fileOut.writelines(lines)
