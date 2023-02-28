PATH_OUTFILES = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\motifs_analysis\\"

class Dataset:
    def __init__(self, motif_file):
        self.__file = motif_file
        self.motifs_list = []
        self.__read_motifs()

    def __read_motifs(self):
        with open(self.__file, "r") as f:
            for line in f:
                motif = line.replace("\n", "").split(",")
                if motif[0] != "X" and motif[2] != "X":
                    self.motifs_list.append(motif)

    def total_motifs(self):
        with open(PATH_OUTFILES + "motifs_analysis_victor2", "w") as f:
            total_paratope_motifs_list = list()
            total_epitope_motifs_list = list()
            motifs_lfr1 = set()
            motifs_lfr2 = set()
            motifs_lfr3 = set()
            motifs_lfr4 = set()
            motifs_cdrl1 = set()
            motifs_cdrl2 = set()
            motifs_cdrl3 = set()
            motifs_hfr1 = set()
            motifs_hfr2 = set()
            motifs_hfr3 = set()
            motifs_hfr4 = set()
            motifs_cdrh1 = set()
            motifs_cdrh2 = set()
            motifs_cdrh3 = set()
            for m in self.motifs_list:
                total_paratope_motifs_list.append(m[0])
                total_epitope_motifs_list.append(m[2])
                if m[1] == "LFR1":
                    motifs_lfr1.add(m[0])
                elif m[1] == "LFR2":
                    motifs_lfr2.add(m[0])
                elif m[1] == "LFR3":
                    motifs_lfr3.add(m[0])
                elif m[1] == "LFR4":
                    motifs_lfr4.add(m[0])
                elif m[1] == "CDR-L1":
                    motifs_cdrl1.add(m[0])
                elif m[1] == "CDR-L2":
                    motifs_cdrl2.add(m[0])
                elif m[1] == "CDR-L3":
                    motifs_cdrl3.add(m[0])
                elif m[1] == "HFR1":
                    motifs_hfr1.add(m[0])
                elif m[1] == "HFR2":
                    motifs_hfr2.add(m[0])
                elif m[1] == "HFR3":
                    motifs_hfr3.add(m[0])
                elif m[1] == "HFR4":
                    motifs_hfr4.add(m[0])
                elif m[1] == "CDR-H1":
                    motifs_cdrh1.add(m[0])
                elif m[1] == "CDR-H2":
                    motifs_cdrh2.add(m[0])
                elif m[1] == "CDR-H3":
                    motifs_cdrh3.add(m[0])

            total_paratope_motifs_set = set(total_paratope_motifs_list)
            total_epitope_motifs_set = set(total_epitope_motifs_list)

            f.write("Motif Analysis\n")
            f.write("Paratope motifs " + str(len(total_paratope_motifs_set)) + "\n")
            f.write("Epitope motifs  " + str(len(total_epitope_motifs_set)) + "\n")
            f.write("CDR - H3 motifs " + str(len(motifs_cdrh3)) + "\n")
            f.write("CDR - H2 motifs " + str(len(motifs_cdrh2)) + "\n")
            f.write("CDR - H1 motifs " + str(len(motifs_cdrh1)) + "\n")

            f.write("CDR - L3 motifs " + str(len(motifs_cdrl3)) + "\n")
            f.write("CDR - L2 motifs " + str(len(motifs_cdrl2)) + "\n")
            f.write("CDR - L1 motifs " + str(len(motifs_cdrl1)) + "\n")

            f.write("HFR 4 motifs    " + str(len(motifs_hfr4)) + "\n")
            f.write("HFR 3 motifs    " + str(len(motifs_hfr3)) + "\n")
            f.write("HFR 2 motifs    " + str(len(motifs_hfr2)) + "\n")
            f.write("HFR 1 motifs    " + str(len(motifs_hfr1)) + "\n")

            f.write("LFR 4 motifs    " + str(len(motifs_lfr4)) + "\n")
            f.write("LFR 3 motifs    " + str(len(motifs_lfr3)) + "\n")
            f.write("LFR 2 motifs    " + str(len(motifs_lfr2)) + "\n")
            f.write("LFR 1 motifs    " + str(len(motifs_lfr1)) + "\n")









ds = Dataset("C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\motifs_akbar\\paratope_motifs_victor2.txt")
ds.total_motifs()

