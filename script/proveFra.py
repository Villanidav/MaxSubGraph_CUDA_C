from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdchem import RWMol
import numpy as np

# largest mapping found so far
incumbent = []
class LabelClass:
    # class used to represent two sets of nodes, from graphs g and h, which share the same label
    # NOTE: the actual value of the label is not needed
    def __init__(self, elems_g, elems_h, rings, adj=0, label=None):
        self.g = elems_g
        self.h = elems_h
        # equals 1 if the label class is adjacent to the current mapping
        self.adj = adj
        self.label = label
        self.rings_g = rings

    def remove(self, graph, elem):
        # remove node and its ring data from label class, from specified graph
        if graph == 0:
            idx = self.g.index(elem)
            del self.g[idx]
            del self.rings_g[idx]
        else:
            idx = self.h.index(elem)
            del self.h[idx]

    # return list with atoms as index and respective ring IDs as data
    def get_ring_match_data(self, elems):

        idxlist = [self.g.index(i) for i in elems]
        res = [self.rings_g[j] for j in idxlist]

        return res


def gen_rotations(s):
    rot_len = len(s)
    tmp = s * 2
    rotations = []

    for i in range(rot_len):
        rotations.append((tmp[i:i + rot_len], i))
    return rotations


# label_classes:    list of all current label classes
def calc_bound(label_classes):
    # maximum possible mapping size = current mapping size +
    # sum of(min(num of nodes with specified label in g, num of nodes with specified label in h)) for each label
    # present in both graphs
    b = 0
    for l in label_classes:
        b += min(len(l.g), len(l.h))
    return b


# label_classes:    list of all current label classes
# map_size:         size of mapping currently being explored
def select_label(label_classes, map_size):
    # selects label form given label classes, based on the maximum number of nodes in either graph with the
    # specified label.
    # if there is no mapping ignore restriction on adjacency, otherwise return only adjacent lable classes since we
    # are looking for connected sub-graphs. Returns none if there are no label classes adjacent to the current mapping

    min_size = 999
    label = None

    for c_label in label_classes:
        if c_label.adj == 1 or map_size == 0:
            c_max_size = max(len(c_label.g), len(c_label.h))
            if c_max_size < min_size:
                min_size = c_max_size
                label = c_label

    return label


# vtx_set:  selected label class
# g:        selected graph
def select_vertex(vtx_set, g):
    # selects node from graph given a label, choosing an adjacent node with the maximum degree

    max_deg = -1
    vtx = 0

    for c_vtx in vtx_set:
        deg = 0

        for i in g[c_vtx]:
            if i != 0:
                deg += 1

        if deg > max_deg:
            max_deg = deg
            vtx = c_vtx

    return vtx


# vtx:  selected node
# g:    selected graph
# edge: bond type
def hood(vtx, g, edge):
    # return the neighbors of a specified node, with the specified bond type.
    friends = []
    for i in range(len(g)):
        if g[i][vtx] == edge and vtx != i:
            friends.append(i)
    return friends



# A ring class is a set of equivalent rings. Two rings are part of the same class if they share the same size and the
# same atom order. Ring equivalence classes are represented here as two arrays, one per molecule,
# containing nothing if the atom is not part of a ring,
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

# l0:   list of all labels in graph 0
# l1:   list of all labels in graph 1
def gen_initial_labels(l0, l1, ring_classes):
    # generates initial label classes based on the set of labels common to both graphs. Initially atoms are
    # grouped only by atom type (if ring match option is enabled, atoms in rings will only be grouped with other atoms
    # in rings of the same type)

    label_classes = []

    common_labels = list(set(l0) & set(l1))

    for label in common_labels:
        # atoms in either molecule which do not have label correspondence are discarded, together with their ring data

        g_elems = [i for i in range(len(l0)) if l0[i] == label]
        g_ring_classes = [ring_classes[a_idx] for a_idx in g_elems]

        h_elems = [j for j in range(len(l1)) if l1[j] == label]

        label_class = LabelClass(g_elems, h_elems, g_ring_classes, label=label)
        label_classes.append(label_class)
    return label_classes

def search_mcs(g0, g1, label_classes, edge_labels, m):
    global incumbent

    # bound calculation
    bound = len(m) + calc_bound(label_classes)

    # if the current mapping is bigger than the biggest one previously found, dethrone the incumbent
    if len(m) > len(incumbent):
        incumbent = m
    # if the incumbent has reached the maximum calculated size given current label classes, return
    if len(incumbent) >= bound:
        return

    # select label form classes
    label_class = select_label(label_classes, len(m))

    # if |m| > 0 and there are no label classes adjacent to m, return
    if label_class is None:
        return

    # select vertex from graph g0
    v = select_vertex(label_class.g, g0)
    # get rings containing v
    v_ring_atoms = label_class.get_ring_match_data([v])[0]

    # cycle through vertices in g1 with selected label
    for idx, w in enumerate(label_class.h):
        # get rings containing w
        # if v and w are not members of at least one shared ring class, v and w cannot be mapped together
        if len(v_ring_atoms) != 0 and (-1 in v_ring_atoms or w not in v_ring_atoms):
            continue
        # label classes for the next recursive call
        l_draft = []
        # create a label class for each type of neighbor (depending on bond and node label) of nodes v and w.
        for label in label_classes:
            for edge_l in edge_labels:
                # v_"conn" and w_"conn" can also contain nodes which are not connected (edge_l == 0)
                v_conn = [vtx for vtx in hood(v, g0, edge_l) if vtx in label.g]
                v_c_rings = label.get_ring_match_data(v_conn)

                w_conn = [vtx for vtx in hood(w, g1, edge_l) if vtx in label.h]

                if len(v_conn) != 0 and len(w_conn) != 0:
                    # the new label class is adjacent to the current mapping if:
                    # 1) its nodes are adjacent to v and w (bond != 0)
                    # 2) its nodes are adjacent to a node pair already present in the mapping
                    adj = (1 if (edge_l != 0 or label.adj == 1) else 0)
                    l_draft.append(LabelClass(v_conn, w_conn, v_c_rings, adj=adj))
        # add (v, w) to current mapping, continue exploring
        search_mcs(g0, g1, l_draft, edge_labels, m + [(v, w)])

    # remove node v from selected label class. If the label class did not contain nodes other than v in graph g0,
    # remove label class
    label_classes.remove(label_class)
    label_class.remove(0, v)

    if len(label_class.g) > 0:
        label_classes.append(label_class)
    # explore consequences of not adding node v to current mapping
    search_mcs(g0, g1, label_classes, edge_labels, m)
def gen_bond_labels(g0, g1):
    # generates bond labels based on the set of bonds common to both graphs.
    common_labels = set([edge for row in g0 for edge in row]) & set([edge for row in g1 for edge in row])
    return list(common_labels)
def mc_split(g0, g1, l0, l1, ring_classes):
    global incumbent
    incumbent = []
    # generate label data
    initial_label_classes = gen_initial_labels(l0, l1, ring_classes)
    edge_labels = gen_bond_labels(g0, g1)

    # search maximum common connected subgraph
    search_mcs(g0, g1, initial_label_classes, edge_labels, [])
    return incumbent
def mol_mcs(mol0, mol1, bond_match=1, ring_match=1, return_map=0):
    mols = [mol0, mol1]

    l0 = [a.GetSymbol() for a in mol0.GetAtoms()]
    l1 = [a.GetSymbol() for a in mol1.GetAtoms()]

    # we need clean labels to reconstruct the molecule in the end
    label_ring_data = [l0.copy(), l1.copy()]

    # if bond_match == 1 the adjacency matrix will contain bond information (i.e: aromatic bond = 1.5, double = 2)
    if bond_match:
        g0 = Chem.rdmolops.GetAdjacencyMatrix(mol0, useBO=True)
        g1 = Chem.rdmolops.GetAdjacencyMatrix(mol1, useBO=True)
    else:
        g0 = Chem.rdmolops.GetAdjacencyMatrix(mol0)
        g1 = Chem.rdmolops.GetAdjacencyMatrix(mol1)

    for line in g1:
        print("\n")
        for idx in line:
            print(idx, end=" ")

    if ring_match:
        # AtomRings() returns a list of rings, represented as a list containing the indexes of atom members
        ring_info = [mol0.GetRingInfo().AtomRings(), mol1.GetRingInfo().AtomRings()]
        for mol_idx in range(2):
            for ring in ring_info[mol_idx]:
                for atom_idx in ring:
                    # if an atom is in 2 rings, make sure not to add a double "R" at the end
                    if label_ring_data[mol_idx][atom_idx][-1] != "R":
                        label_ring_data[mol_idx][atom_idx] += "R"
                     

    ring_classes = gen_ring_classes(mol0, mol1)

    mapping = mc_split(g0, g1, label_ring_data[0], label_ring_data[1], ring_classes)
    return mapping
    
def smiles_mcs(s0, s1, bond_match=1, ring_match=1):

    mol0 = Chem.MolFromSmiles(s0)
    mol1 = Chem.MolFromSmiles(s1)

    return mol_mcs(mol0, mol1, bond_match, ring_match)

"""

l0 = ["C", "C", "O", "C", "H", "H", "O"]
l1 = ["C","H","F","O","H","C","B","R"]

ring_info_m0 = [ [0, 1, 2], [4, 3, 5]]
ring_info_m1  = [[4, 1, 0],  [1, 2, 3], [1, 4, 5], [5, 0, 3]]


gen_classes = gen_ring_classes(l0,l1,ring_info_m0,ring_info_m1)
i = 0

#for posizione in gen_classes:
#    print("\nidx:", i)
#    i += 1

#    for j in posizione:
#        print(j, end=" ")

initial_label = gen_initial_labels(l0,l1,gen_classes)

#
#    for label in initial_label:
#        print("\n" + label.label + ":")
#        print("\ng:")
#        for idx in label.g :
#            print(idx)
#        print("\nh:")
#        for idx in label.h:
#            print(idx)
#
#
"""

smile0 = "CN(c1ccc(cc1)c2nnn(CC(=O)Nc3ccc4nc(oc4c3)c5ccccc5Cl)n2)c6cc(C)c(N)cn6"
smile1 = "O=C(Cn1nnc(n1)c2ccc(Nc3ccccn3)cc2)Nc4ccc5nc(oc5c4)c6ccccc6"
result = smiles_mcs(smile0,smile1,bond_match=1, ring_match=1)

print(result)

