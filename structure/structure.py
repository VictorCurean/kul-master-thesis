import Bio.PDB.Selection
from utils import *
from Bio.PDB import NeighborSearch, PDBParser
import re

parser = PDBParser(PERMISSIVE=False)

def get_structure(pdb_id, file):
    return parser.get_structure(id=pdb_id, file=file)

def get_opposite_complex(origin, targets):
    type_origin = origin.get_full_id()[2]
    targets_to_keep = list()
    for t in targets:
        #origin is antigen, target in heavy- or light chain:
        if t.get_full_id()[2] in ["H", "L"] and type_origin == "C":
            targets_to_keep.append(t)
        #origin in heavy- or light chain, target in antigen
        elif t.get_full_id()[2] == "C" and type_origin in ["H", "L"]:
            targets_to_keep.append(t)

    return set(targets_to_keep)

def get_neighboring_residues(struct, target_residue, cutoff):
    atoms = list(struct.get_atoms())
    pdb_nbsearch = Bio.PDB.NeighborSearch(atoms)
    res_atoms = target_residue.get_atoms()

    neighboring_residues = set()
    for a in list(res_atoms):
        neighboring_residues |= set(pdb_nbsearch.search(a.get_coord(), cutoff, level='R'))

    return get_opposite_complex(target_residue, neighboring_residues)

def get_residues_by_type(struct, type):
    # H - heavy chain, L - light chain, C - antigen
    return [r for r in struct.get_residues() if r.get_full_id()[2] == type]

def write_annotated_sequence(file, struct, type):
    SEQ = ""
    ANN = ""
    for r in get_residues_by_type(struct, type):
        SEQ += AA_3to1[r.get_resname()]
        nb_res = get_neighboring_residues(struct, r, cutoff=5)
        if len(nb_res) != 0:
            ANN += "#"
        else:
            ANN += "_"
    file.write(SEQ + "\n")
    file.write(ANN + "\n")

    if type in ["H", "L"]:
        file.write(annotate_sequence_by_region(SEQ, type) + "\n")

def annotate_sequence_by_region(seq, type):
    ann = ""
    res = get_martin_numbering(seq)
    for i in range(len(seq)):
        pos_martin = int(re.sub("[^0-9]", "", res[i].split(" ")[0]))
        if type == "H":
            if pos_martin in MARTIN_NUMBERING_HFR:
                ann += "_"
            else:
                ann += "C"

        elif type == "L":
            if pos_martin in MARTIN_NUMBERING_LFR:
                ann += "_"
            else:
                ann += "C"
    return ann

def annotate_sequence(struct, pdb_id):
    OUT_DIR = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\annotated_sequences\\"
    with open(OUT_DIR+pdb_id+".seqann", "w") as file:
        file.write("ID\t" + pdb_id + "\n")

        #LIGHT CHAIN
        file.write("Light chain" + "\n")
        write_annotated_sequence(file, struct, "L")

        #HEAVY CHAIN
        file.write("Heavy chain" + "\n")
        write_annotated_sequence(file, struct, "H")

        #ANTIGEN
        file.write("Antigen" + "\n")
        write_annotated_sequence(file, struct, "C")

def get_motifs():
    pass

struct = get_structure("1A2Y", "C:\\Users\\curea\\PycharmProjects\\kul-thesis-ab\\datasets\\1A2Y_1.pdb")

for r in struct.get_residues():
    a = get_neighboring_residues(struct, r, cutoff=5)
    if len(a) != 0:
        b = 1


annotate_sequence(struct, "1A2Y")

