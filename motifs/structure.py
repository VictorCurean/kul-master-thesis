import Bio.PDB.Selection
from utils import *
from Bio.PDB import NeighborSearch, PDBParser
import re
from os import listdir
from os.path import isfile, join

parser = PDBParser(PERMISSIVE=False)

def get_structure(pdb_id, file):
    return parser.get_structure(id=pdb_id, file=file)

def get_opposite_complex(origin, targets):
    type_origin = origin.get_full_id()[2]
    targets_to_keep = list()
    for t in targets:
        #origin is antigen, target in heavy- or light chain:
        if t.get_full_id()[2] in ["H", "L"] and type_origin not in ["H", "L"]:
            targets_to_keep.append(t)
        #origin in heavy- or light chain, target in antigen
        elif t.get_full_id()[2] not in ["H", "L"] and type_origin in ["H", "L"]:
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

    if type in ["H", "L"]:
        return [r for r in struct.get_residues() if r.get_full_id()[2] == type]
    else:
        return [r for r in struct.get_residues() if r.get_full_id()[2] not in ["H", "L"]]

def translate_to_1_notation(seq):
    res = ""
    for r in seq:
        res += AA_3to1[r.get_resname()]
    return res

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

    return SEQ, ANN

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

# struct = get_structure("1A2Y", "C:\\Users\\curea\\PycharmProjects\\kul-thesis-ab\\datasets\\1A2Y_1.pdb")
#
# for r in struct.get_residues():
#     a = get_neighboring_residues(struct, r, cutoff=5)
#     if len(a) != 0:
#         b = 1


# annotate_sequence(struct, "1A2Y")

class Complex:
    def __init__(self, pdb_id, pdb_file):
        self.__pdb_id = pdb_id
        self.__struct = get_structure(pdb_id, pdb_file)
        self.__lc = translate_to_1_notation(get_residues_by_type(self.__struct, "L"))
        self.__hc = translate_to_1_notation(get_residues_by_type(self.__struct, "H"))
        self.martin_lc = get_martin_numbering(self.__lc)
        self.martin_hc = get_martin_numbering(self.__hc)
        self.__cutoff = 5
        self.motifs = []

    def get_struct(self):
        return self.__struct

    def __get_martin_scheme_by_type(self, type):
        if type == "L":
            return [x for x in self.martin_lc if x != '']
        elif type == "H":
            return [x for x in self.martin_hc if x != '']


    def __get_motifs_para_epi(self, begin, end, type, loc):
        ab_residues = get_residues_by_type(self.__struct, type)
        antigen_residues = get_residues_by_type(self.__struct, "C")

        antigen_nb_residues = []
        ab_nb_residues = []

        martin_scheme = self.__get_martin_scheme_by_type(type)

        for i in range(len(ab_residues)):
            martin_pos = int(re.sub("[^0-9]", "", martin_scheme[i].split(" ")[0]))
            if begin <= martin_pos <= end:
                nb = get_neighboring_residues(self.__struct, ab_residues[i], self.__cutoff)
                if len(nb) != 0:
                    ab_nb_residues.append(ab_residues[i])
                    antigen_nb_residues += nb


        antigen_nb_residues = list(set(antigen_nb_residues))
        antigen_nb_residues.sort(key=lambda x : x.id[1], reverse=False)

        if len(antigen_nb_residues + ab_nb_residues) != 0:
            ant_res, ant_enc = self.__get_gaps(antigen_nb_residues, antigen_residues)
            lc_res, lc_enc = self.__get_gaps(ab_nb_residues, ab_residues)

            self.motifs.append(Motif(self.__pdb_id, "para-epi", lc_res, ant_res, lc_enc, ant_enc, loc))


    def __get_gaps(self, res_list, full_res):
        residues_no_gaps = []
        bool_encodings = []

        min_res = min(res_list, key=lambda res : res.id[1])
        max_res = max(res_list, key=lambda res: res.id[1])

        index_list = [r.id[1] for r in full_res if min_res.id[1] <= r.id[1] <= max_res.id[1]]

        for i in index_list:
            if i in [x.id[1] for x in res_list]:
                residues_no_gaps.append([x for x in res_list if x.id[1] == i][0])
                bool_encodings.append(True)
            else:
                residues_no_gaps.append([x for x in full_res if x.id[1] == i][0])
                bool_encodings.append(False)

        return residues_no_gaps, bool_encodings

    def get_motifs(self):
        #LFR1
        self.__get_motifs_para_epi(1, 23, "L", "LFR1")
        #CDR-L1
        self.__get_motifs_para_epi(24, 34, "L", "CDR-L1")
        #LFR2
        self.__get_motifs_para_epi(35, 49, "L", "LFR2")
        #CDR-L2
        self.__get_motifs_para_epi(50, 56, "L", "CDR-L2")
        #LFR3
        self.__get_motifs_para_epi(57, 88, "L", "LFR3")
        #CDR-L3
        self.__get_motifs_para_epi(89, 97, "L", "CDR-L3")
        #LFR4
        self.__get_motifs_para_epi(98, 110, "L", "LFR4")

        #HFR1
        self.__get_motifs_para_epi(1, 30, "H", "HFR1")
        #CDR-H1
        self.__get_motifs_para_epi(31, 35, "H", "CDR-H1")
        #HFR2
        self.__get_motifs_para_epi(36, 49, "H", "HFR2")
        #CDR-H2
        self.__get_motifs_para_epi(50, 65, "H", "CDR-H2")
        #HFR3
        self.__get_motifs_para_epi(66, 94, "H", "HFR3")
        #CDR-H3
        self.__get_motifs_para_epi(95, 102, "H", "CDR-H3")
        #HFR4
        self.__get_motifs_para_epi(103, 113, "H", "HFR4")

    def __cleanup_motifs(self):
        self.motifs = [m for m in self.motifs if m.is_trivial() is False]

    def write_motifs(self, paratope_file):
        for m in self.motifs:
            m.write_motifs_akbar_encoding(paratope_file)

    def write_motifs_victor1(self, paratope_file):
        for m in self.motifs:
            m.write_motifs_victor1_encoding(paratope_file)

    def write_motifs_victor2(self, paratope_file):
        for m in self.motifs:
            m.write_motifs_victor2_encoding(paratope_file)

    def check_motifs(self):
        if len(self.motifs) == 0:
            print(self.__pdb_id + ".... no motifs")


class Motif:
    def __init__(self, pdb_id, type, res_origin, res_target, encodings_origin, encodings_target, ab_location):
        self.pdb_id = pdb_id
        self.type = type
        self.res_origin = res_origin
        self.res_target = res_target
        self.encodings_origin = encodings_origin
        self.encodings_target = encodings_target
        self.ab_location = ab_location

    def write_motifs_akbar_encoding(self, paratope_file):
        notation_origin = self.__get_akbar_notation_from_res(self.encodings_origin)
        notation_target = self.__get_akbar_notation_from_res(self.encodings_target)

        with open(paratope_file, "a") as pf:
            pf.write(notation_origin + "," + self.ab_location + "," + notation_target + "," + self.pdb_id + "\n")


    def write_motifs_victor1_encoding(self, paratope_file):
        notation_origin = self.__get_victor1_encoding_from_res(self.encodings_origin, self.res_origin)
        notation_target = self.__get_victor1_encoding_from_res(self.encodings_target, self.res_target)

        with open(paratope_file, "a") as pf:
            pf.write(notation_origin + "," + self.ab_location + "," + notation_target + "\n")

    def write_motifs_victor2_encoding(self, paratope_file):
        notation_origin = self.__get_victor2_encoding_from_res(self.encodings_origin, self.res_origin)
        notation_target = self.__get_victor2_encoding_from_res(self.encodings_target, self.res_target)

        with open(paratope_file, "a") as pf:
            pf.write(notation_origin + "," + self.ab_location + "," + notation_target + "\n")
    def __get_akbar_notation_from_res(self, res_encodings):
        notation = ""
        counter = 0
        for e in res_encodings:
            if e is True:
                if counter != 0:
                    notation += str(counter)
                    counter = 0
                notation += "X"
            else:
                counter += 1

        return notation

    def __get_victor1_encoding_from_res(self, res_encodings, res_names):
        notation = ""
        counter = 0
        for i in range(len(res_encodings)):
            if res_encodings[i] is True:
                if counter != 0:
                    notation += str(counter)
                    counter = 0
                notation += AA_TO_TYPE[AA_3to1[res_names[i].resname]]
            else:
                counter += 1

        return notation

    def __get_victor2_encoding_from_res(self, res_encodings, res_names):
        notation = ""
        counter = 0
        for i in range(len(res_encodings)):
            if res_encodings[i] is True:
                if counter != 0:
                    notation += str(counter)
                    counter = 0
                notation += AA_TO_TYPE[AA_3to1[res_names[i].resname]]
            else:
                counter += 1

        return notation

    def is_trivial(self):
        if len(self.res_origin) == 1:
            return True
        return False


# a = Complex("1HI6", "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\datasets\\1HI6_1.pdb")
# a.get_motifs()

paratope_file = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\motifs_akbar\\motifs_akbar.txt"
#paratope_file2 = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\motifs_akbar\\paratope_motifs_victor2.txt"
# a.write_motifs(paratope_file, epitope_file)

def read_all_motifs():
    path = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\datasets"

    open(paratope_file, "w").close()
    files = [f for f in listdir(path) if isfile(join(path, f))]

    err_files = 0

    for f in files:
        print(f)
        comp = Complex(f, path + "\\" + f)
        try:
            comp.get_motifs()
            comp.check_motifs()
            comp.write_motifs(paratope_file)
        except IndexError:
            print(f + " .... Index Error")
            err_files += 1
    print(err_files)

read_all_motifs()
