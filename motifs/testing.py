from structure import get_structure, get_neighboring_residues, annotate_sequence, Complex

# struct = get_structure("1A14", "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\datasets\\1A14_1.pdb")
# res = list(struct.get_residues())
# chains = list(struct.get_chains())
# for r in struct.get_residues():
#     a = get_neighboring_residues(struct, r, cutoff=5)
#     if len(a) != 0:
#         b = 1
#
#
# annotate_sequence(struct, "1A14")

path = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\datasets"
f = "1FNS_1.pdb"

comp = Complex(f, path + "\\" + f)
comp.get_motifs()