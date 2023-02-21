from structure import get_structure, get_neighboring_residues, annotate_sequence


struct = get_structure("1A14", "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\datasets\\1A14_1.pdb")
res = list(struct.get_residues())
chains = list(struct.get_chains())
for r in struct.get_residues():
    a = get_neighboring_residues(struct, r, cutoff=5)
    if len(a) != 0:
        b = 1


annotate_sequence(struct, "1A14")