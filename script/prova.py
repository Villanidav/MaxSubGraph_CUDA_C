from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdchem import RWMol
import numpy as np

def gen_rotations(s):
    rot_len = len(s)
    tmp = s * 2
    rotations = []

    for i in range(rot_len):
        rotations.append((tmp[i:i + rot_len], i))
    return rotations


def gen_ring_classes_s(l0, l1, ring_info_m0, ring_info_m1 ):
    ring_comp_m0 = [[] for _ in l0]
    
    for r0 in ring_info_m0:
        for a0 in r0:
            ring_comp_m0[a0] = [-1]

    # iterate through rings to match them to their compatibility class, marking each atom as a member
    for r0 in ring_info_m0:

        r0_label = ""

        for a0 in r0:
            r0_label += str(l0[a0])
        r0_label_rev = r0_label[::-1]

        for r1 in ring_info_m1:

            if len(r0) == len(r1):
                r1_label = ""

                for a1 in r1:
                    r1_label += str(l1[a1])

                # to check if the relative order of atoms is the same between two rings, check if the strings
                # composed by the member atoms are rotations of each other. Since only symmetrical rings
                # can have "equivalent" atoms, we don't need to reverse the ring label

                r0_rots = [rot for rot in gen_rotations(r0_label) if rot[0] == r1_label]
                print(r0_rots)
                r0_rev_rots = [rot for rot in gen_rotations(r0_label_rev) if rot[0] == r1_label]
                print(r0_rev_rots)
                i = 0 
                for (label, amt) in r0_rots:
                    for idx, a0 in enumerate(r0):
                        print(r0[idx])
                        if ring_comp_m0[a0] == [-1]:
                            ring_comp_m0[a0] = [r1[(idx - amt)]]
                            print(i ,"la posizione idx dell'anello r1 considerata:", idx - amt, "ed è ", r1[idx-amt])
                            print(ring_comp_m0)
                            i = i+1
                        else:
                            ring_comp_m0[a0].append(r1[(idx - amt)])
                            print(i ,"la posizione idx dell'anello r1 considerata:", idx - amt, "ed è ", r1[idx-amt])
                            print(ring_comp_m0)
                            i = i+1

                for (label, amt) in r0_rev_rots:
                    for idx, a0 in enumerate(r0):
                        print(r0[idx])
                        rev_idx = len(r1) - idx - 1
                        if ring_comp_m0[a0] == [-1]:
                            ring_comp_m0[a0] = [r1[(rev_idx - amt)]]
                            print(i , "la posizione idx dell'anello r1 considerata: " , rev_idx - amt, "ed è : ", r1[rev_idx - amt])
                            print(ring_comp_m0)
                            i = i+1
                        elif r1[(rev_idx - amt)] not in ring_comp_m0[a0]:
                            ring_comp_m0[a0].append(r1[(rev_idx - amt)])
                            print(i , "la posizione idx dell'anello r1 considerata: " , rev_idx - amt, "ed è : ", r1[rev_idx - amt])
                            print(ring_comp_m0)
                            i = i+1
    return ring_comp_m0

def s_smiles_mcs(s0, s1):

    mol0 = Chem.MolFromSmiles(s0)
    mol1 = Chem.MolFromSmiles(s1)

    return gen_ring_classes(mol0,mol1)

def gen_ring_classes(mol0, mol1):




    l0 = [a.GetSymbol() for a in mol0.GetAtoms()]
    l1 = [a.GetSymbol() for a in mol1.GetAtoms()]

    # arrays containing ring class information for each of the molecule's atoms
    ring_comp_m0 = [[] for _ in l0]

    # AtomRings() returns a list of rings from the molecule containing the indexes of atom members
    ring_info_m0 = mol0.GetRingInfo().AtomRings()
    ring_info_m1 = mol1.GetRingInfo().AtomRings()

    for r0 in ring_info_m0:
        for a0 in r0:
            ring_comp_m0[a0] = [-1]

    # iterate through rings to match them to their compatibility class, marking each atom as a member
    for r0 in ring_info_m0:

        r0_label = ""

        for a0 in r0:
            r0_label += str(l0[a0])
        r0_label_rev = r0_label[::-1]

        for r1 in ring_info_m1:

            if len(r0) == len(r1):
                r1_label = ""

                for a1 in r1:
                    r1_label += str(l1[a1])

                # to check if the relative order of atoms is the same between two rings, check if the strings
                # composed by the member atoms are rotations of each other. Since only symmetrical rings
                # can have "equivalent" atoms, we don't need to reverse the ring label

                r0_rots = [rot for rot in gen_rotations(r0_label) if rot[0] == r1_label]
                r0_rev_rots = [rot for rot in gen_rotations(r0_label_rev) if rot[0] == r1_label]

                for (label, amt) in r0_rots:
                    for idx, a0 in enumerate(r0):
                        if ring_comp_m0[a0] == [-1]:
                            ring_comp_m0[a0] = [r1[(idx - amt)]]
                        else:
                            ring_comp_m0[a0].append(r1[(idx - amt)])

                for (label, amt) in r0_rev_rots:
                    for idx, a0 in enumerate(r0):
                        rev_idx = len(r1) - idx - 1
                        if ring_comp_m0[a0] == [-1]:
                            ring_comp_m0[a0] = [r1[(rev_idx - amt)]]
                        elif r1[(rev_idx - amt)] not in ring_comp_m0[a0]:
                            ring_comp_m0[a0].append(r1[(rev_idx - amt)])
    return ring_comp_m0


'''
l0 = ["C","C","O","C","H","H","O"]
l1 = ["C","H","F","O","H","C","B","R"]
ring_info_m0 = [[0,1,2],[4,3,5]]
ring_info_m1 = [[4,1,0],[1,2,3],[1,4,5],[5,0,3]]


generated_rings = gen_ring_classes_s(l0,l1,ring_info_m0,ring_info_m1)

pos = 0
for vect in generated_rings:
    print("\n: posizione:", pos, "\n")
    pos += 1
    for i in vect:
        print(i, end=" ")
'''

smile0 = "CN(c1ccc(cc1)c2nnn(CC(=O)Nc3ccc4nc(oc4c3)c5ccccc5Cl)n2)c6cc(C)c(N)cn6"
smile1 = "O=C(Cn1nnc(n1)c2ccc(Nc3ccccn3)cc2)Nc4ccc5nc(oc5c4)c6ccccc6"


generated_rings = s_smiles_mcs(smile0,smile1)

pos = 0
for vect in generated_rings:
    print("\n: posizione:", pos, "\n")
    pos += 1
    for i in vect:
        print(i, end=" ")






