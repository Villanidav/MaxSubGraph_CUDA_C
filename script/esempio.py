
filename = "smiles.txt"
smiles = []

try:
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip().split()
            first = line[1]  # First string after skipping the first int
            second = line[5]  # Second string after skipping four ints
            smiles.append((first, second))
except FileNotFoundError:
    print("Error opening file:", filename)



# No need to explicitly close the file in Python, as it's handled by the "with" statement

# Open the file in write mode (creates the file if it doesn't exist)
with open("outputPython.txt", "w") as f:
    for i in range(len(smiles) -1):     
        result = smiles_mcs(smiles[i][0],smiles[i][1],bond_match=1, ring_match=1)
        l0 = [a.GetSymbol() for a in result.GetAtoms()]
        print(l0, file=f)
            
    







