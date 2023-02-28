from os import listdir
from os.path import isfile, join

PATH = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\"
path_pdb_files = PATH + "datasets\\"
path_motifs = PATH + "motifs_akbar\\"
path_complex_analysis = PATH + "foldx\\out.complex_analysis"


def build_datamatrix():
    files = [f for f in listdir(path_pdb_files) if isfile(join(path_pdb_files, f))]
    for f in files:
        print(f)

def get_nr_connecting_residues(pdb_id):
    pass

def get_ddg_complex_analysis(pdb_id):
    pass

def get_no_of_motifs_epitope(pdb_id):
    pass

def get_no_of_motifs_paratope(pdb_id):
    pass

def get_no_of_continous_motifs_epitope(pdb):
    pass

def get_no_of_continous_motifs_paratope(pdb):
    pass

build_datamatrix()