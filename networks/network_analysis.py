import networkx
import networkx as nx


PATH = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\motifs_akbar\\"
def get_data_from_file(file):
    paratope_nodes = dict()
    epitope_nodes = dict()
    edges = dict()

    with open(file, "r") as f:
        for line in f:
            motif = line.replace("\n", "").split(",")
            para_motif = "P" + motif[0]
            epi_motif = "E" + motif[2]

            if para_motif in paratope_nodes.keys():
                paratope_nodes[para_motif] += 1
            elif motif[0] not in paratope_nodes.keys():
                paratope_nodes[para_motif] = 1

            if epi_motif in epitope_nodes.keys():
                epitope_nodes[epi_motif] += 1
            elif epi_motif not in epitope_nodes.keys():
                epitope_nodes[epi_motif] = 1

            if (para_motif, epi_motif,) in edges.keys():
                edges[(para_motif, epi_motif,)] += 1
            elif (para_motif, epi_motif,) not in edges.keys():
                edges[(para_motif, epi_motif,)] = 1

    return paratope_nodes, epitope_nodes, edges

def build_network(paratope_nodes, epitope_nodes, edges):
    reactivity_network = nx.Graph()
    for pn in paratope_nodes.keys():
        reactivity_network.add_node(pn, connections=paratope_nodes[pn], bipartite=0, type="paratope")
    for en in epitope_nodes.keys():
        reactivity_network.add_node(en, connections=epitope_nodes[en], bipartitie=1, type="epitope")
    for e in edges.keys():
        reactivity_network.add_edge(e[0], e[1], weight=edges[e])

    return reactivity_network


paratope_nodes, epitope_nodes, edges = get_data_from_file(PATH + "motifs_final.txt")
network = build_network(paratope_nodes, epitope_nodes, edges)
network.remove_node("PX")
network.remove_node("EX")

#networkx.write_graphml(network, PATH+"motifs_cytoscape.xml")

nodes = list(nx.nodes(network))
degree_hist = sorted(nx.degree_histogram(network), reverse=True)
b = 1