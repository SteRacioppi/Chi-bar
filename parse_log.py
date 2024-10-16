import re
import argparse

# Function to parse alpha electrons
def parse_alpha_occupations(lines):
    alpha_occup = []
    float_regex = re.compile(r'[-+]?\d*\.\d+|\d+')
    
    alpha = False
    for line in lines:
        if "Alpha  occ. eigenvalues" in line:
            alpha = True
            eigenvalues = [float(x) for x in float_regex.findall(line)]
            alpha_occup.extend(eigenvalues)
        elif "Alpha virt. eigenvalues" in line:
            # Stop parsing at the start of virtual orbitals
            break
        elif alpha:
            # Additional occupied orbitals continued on the next lines
            eigenvalues = [float(x) for x in float_regex.findall(line)]
            alpha_occup.extend(eigenvalues)
    
    return alpha_occup

# Function to parse beta electrons
def parse_beta_occupations(lines):
    beta_occup = []
    float_regex = re.compile(r'[-+]?\d*\.\d+|\d+')
    
    beta = False
    for line in lines:
        if "Beta  occ. eigenvalues" in line:
            beta = True
            eigenvalues = [float(x) for x in float_regex.findall(line)]
            beta_occup.extend(eigenvalues)
        elif "Beta virt. eigenvalues" in line:
            # Stop parsing at the start of virtual orbitals
            break
        elif beta:
            # Additional occupied orbitals continued on the next lines
            eigenvalues = [float(x) for x in float_regex.findall(line)]
            beta_occup.extend(eigenvalues)

    return beta_occup

# Function to print the orbitals
def print_occupied_orbitals(alpha_occup, beta_occup):
    # Print Alpha orbitals
    print("Alpha Occupied Molecular Orbitals (with index):")
    for i, orb in enumerate(alpha_occup, 1):
        print(f"Orbital {i}: {orb:.6f}")

    # Print Beta orbitals
    if beta_occup:
        print("\nBeta Occupied Molecular Orbitals (with index):")
        for i, orb in enumerate(beta_occup, 1):
            print(f"Orbital {i}: {orb:.6f}")
    else:
        print("\nNo Beta occupied orbitals found.")

# Main function to handle parsing and output
def main():
    parser = argparse.ArgumentParser(description="Parse a Gaussian log file and print the occupied molecular orbitals.")
    parser.add_argument('logfile', type=str, help='Path to the Gaussian log file')
    
    args = parser.parse_args()

    # Read the file
    with open(args.logfile, 'r') as f:
        lines = f.readlines()

    # Parse alpha and beta occupations separately
    alpha_occup = parse_alpha_occupations(lines)
    beta_occup = parse_beta_occupations(lines)

    # Sort in ascending order (core to valence)
    alpha_occup.sort()
    beta_occup.sort()

    # Print the results
    print_occupied_orbitals(alpha_occup, beta_occup)

if __name__ == "__main__":
    main()


