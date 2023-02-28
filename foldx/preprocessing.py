from os import listdir
from os.path import isfile, join

PATH = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\"
path_pdb_files = PATH + "datasets\\"
path_motifs = PATH + "motifs_akbar\\"
path_complex_analysis = PATH + "foldx\\out.complex_analysis\\"
path_datamatrix = PATH + "foldx\\datamatrix\\"


def build_datamatrix():
    files = [f for f in listdir(path_pdb_files) if isfile(join(path_pdb_files, f))]
    outfile = path_datamatrix + "datamatrix.txt"

    with open(outfile, "w") as of:
        of.write("pdb\tno_motifs\tconn_res_paratope\tconn_res_epitope\tenergy\thfr1\thfr2\thfr3\thfr4\tlfr1\tlfr2\tlfr3\tlfr4\tcdrh1\tcdrh2\tcdrh3\tcdrl1\tcdrl2\tcdrl3\tavg_len_para\tavg_len_epi\n")
        for f in files:
            pdb_id = f.split(".")[0]
            no_of_motifs, conn_res_paratope, conn_res_epitope, \
            motifs_hfr1, motifs_hfr2, motifs_hfr3, motifs_hfr4, \
            motifs_lfr1, motifs_lfr2, motifs_lfr3, motifs_lfr4, \
            motifs_cdrh1, motifs_cdrh2, motifs_cdrh3, \
            motifs_cdrl1, motifs_cdrl2, motifs_cdrl3, \
            avg_length_motif_epitope, avg_length_motif_paratope = get_nr_connecting_residues(pdb_id)
            energy = get_ddg_complex_analysis(pdb_id)

            of.write(pdb_id + "\t" + str(no_of_motifs) + "\t" + str(conn_res_paratope) + "\t" +
                     str(conn_res_epitope) + "\t" + str(energy) + "\t" + str(motifs_hfr1) + "\t" + str(motifs_hfr2) +
                     "\t" + str(motifs_hfr3) + "\t" + str(motifs_hfr4) + "\t" +
                     str(motifs_lfr1) + "\t" + str(motifs_lfr2) + "\t" + str(motifs_lfr3) + "\t" + str(motifs_lfr4) + "\t" +
                     str(motifs_cdrh1) + "\t" + str(motifs_cdrh2) + "\t" + str(motifs_cdrh3) + "\t" +
                     str(motifs_cdrl1) + "\t" + str(motifs_cdrl2) + "\t" + str(motifs_cdrl3) + "\t" +
                     str(avg_length_motif_paratope) + "\t" + str(avg_length_motif_epitope) + "\n")



def get_nr_connecting_residues(pdb_id):
    motifs = list()
    with open(path_motifs + "motifs_final_pdb.txt", "r") as f:
        for line in f:
            motif = line.split(",")
            pdb = motif[3].strip("\n").split(".")[0]
            if pdb == pdb_id:
                motifs.append(motif)

    conn_res_paratope = 0
    conn_res_epitope = 0
    no_of_motifs = 0
    motifs_cdrh1 = 0
    motifs_cdrh2 = 0
    motifs_cdrh3 = 0

    motifs_cdrl1 = 0
    motifs_cdrl2 = 0
    motifs_cdrl3 = 0

    motifs_hfr1 = 0
    motifs_hfr2 = 0
    motifs_hfr3 = 0
    motifs_hfr4 = 0

    motifs_lfr1 = 0
    motifs_lfr2 = 0
    motifs_lfr3 = 0
    motifs_lfr4 = 0


    for m in motifs:
        no_of_motifs += 1
        conn_res_paratope += m[0].count("X")
        conn_res_epitope += m[2].count("X")

        if m[1] == "CDR-H1":
            motifs_cdrh1 += 1
        elif m[1] == "CDR-H2":
            motifs_cdrh2 += 1
        elif m[1] == "CDR-H3":
            motifs_cdrh3 += 1

        elif m[1] == "CDR-L1":
            motifs_cdrl1 += 1
        elif m[1] == "CDR-L2":
            motifs_cdrl2 += 1
        elif m[1] == "CDR-L3":
            motifs_cdrl3 += 1

        elif m[1] == "HFR1":
            motifs_hfr1 += 1
        elif m[1] == "HFR2":
            motifs_hfr2 += 1
        elif m[1] == "HFR3":
            motifs_hfr3 += 1
        elif m[1] == "HFR4":
            motifs_hfr4 += 1

        elif m[1] == "LFR1":
            motifs_lfr1 += 1
        elif m[1] == "LFR2":
            motifs_lfr2 += 1
        elif m[1] == "LFR3":
            motifs_lfr3 += 1
        elif m[1] == "LFR4":
            motifs_lfr4 += 1

        avg_length_motif_paratope = conn_res_paratope/no_of_motifs
        avg_length_motif_epitope = conn_res_epitope/no_of_motifs

    return no_of_motifs, conn_res_paratope, conn_res_epitope, \
           motifs_hfr1, motifs_hfr2, motifs_hfr3, motifs_hfr4, \
           motifs_lfr1, motifs_lfr2, motifs_lfr3, motifs_lfr4, \
           motifs_cdrh1, motifs_cdrh2, motifs_cdrh3,\
           motifs_cdrl1, motifs_cdrl2, motifs_cdrl3,\
           avg_length_motif_epitope, avg_length_motif_paratope



def get_ddg_complex_analysis(pdb_id):
    file = path_complex_analysis + "Summary_" + pdb_id + "_AC.fxout"
    with open(file, "r") as f:
        for i, line in enumerate(f):
            if i == 9:
                return line.split("\t")[5]

build_datamatrix()

#get_ddg_complex_analysis("2JEL_1")