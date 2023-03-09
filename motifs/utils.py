import requests


AA_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
     'H1S': 'H', 'H2S': 'H',
    }


MARTIN_NUMBERING = {"LFR1": (1, 23), "CDR-L1": (24, 34), "LFR2": (35, 49), "CDR-L2": (50, 56), "LFR3": (57, 88),
                    "CDR-L3": (89, 97), "LFR4": (98, 110), "HFR1": (1, 30), "CDR-H1": (31, 35), "HFR2": (36, 49),
                    "CDR-H2": (50, 65), "HFR3": (66, 94), "CDR-H3": (95, 102), "HFR4": (103, 113)}

MARTIN_NUMBERING_LFR = list(range(1, 24)) + list(range(35, 50)) + list(range(57, 89)) + list(range(98, 111))
MARTIN_NUMBERING_HFR = list(range(1, 31)) + list(range(36, 50)) + list(range(66, 95)) + list(range(103, 113))


AA_TO_TYPE = {"G": "A", "A": "A", "V": "A", "L": "A", "I": "A", "P": "A", "M": "A",
              "F": "B", "Y": "B", "W": "B",
              "S": "C", "T": "C", "C": "C", "N": "C", "Q": "C",
              "D": "D", "E": "D",
              "R": "E", "H": "E", "K": "E"}


AA_TO_TYPE_LOOSE = {"G": "Z", "A": "Z", "V": "Z", "L": "Z", "I": "Z", "P": "Z", "M": "Z",
              "F": "Z", "Y": "Z", "W": "Z",
              "S": "Y", "T": "Y", "C": "Y", "N": "Y", "Q": "Y",
              "D": "Y", "E": "Y",
              "R": "Y", "H": "Y", "K": "Y"}

def get_martin_numbering(seq):
    API_NUMERBING_URL = "http://www.bioinf.org.uk/abs/abnum/abnum.cgi?plain=1&aaseq={seq}&scheme={scheme}"

    response = requests.get(API_NUMERBING_URL.format(seq=seq, scheme="-m"))
    return response.text.split("\n")

#get_martin_numbering("QVQLQESGPGLVAPSQSLSITCTVSGFSLTGYGVNWVRQPPGKGLEWLGMIWGDGNTDYNSALKSRLSISKDNSKSQVFLKMNSLHTDDTARYYCARERDYRLDYWGQGTTLTVSS")