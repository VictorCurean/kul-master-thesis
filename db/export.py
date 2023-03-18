import itertools

import psycopg2
import csv

OUT_PATH = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\db\\csv\\"


DB_NAME = "kul_thesis"
USER = "postgres"
PASSWORD = "admin"
HOST = "localhost"

def get_conn():
    return psycopg2.connect(host=HOST, database=DB_NAME, user=USER, password=PASSWORD)

def get_interaction_values():
    conn = get_conn()
    cur = conn.cursor()

    SELECT_QUERY = "SELECT * FROM complex_analysis_interaction"
    cur.execute(SELECT_QUERY)
    result = cur.fetchall()
    with open(OUT_PATH + "interaction_data.csv", "w") as f:
        out_file = csv.writer(f, lineterminator="\n")
        column_names = [i[0] for i in cur.description]
        out_file.writerow(column_names)
        out_file.writerows(result)


def get_energy_distribution_for_motif(para_motif_heavy_chain, location):
    conn = get_conn()
    cur = conn.cursor()

    SELECT_QUERY = "SELECT * FROM motifs WHERE para_motif LIKE '{para_motif}' AND location LIKE '{location}'"

    cur.execute(SELECT_QUERY.format(para_motif=para_motif_heavy_chain, location=location))
    res = cur.fetchall()

    results = list()
    for a in res:

        QUERY = """
        SELECT ala_scan_energy 
        FROM residues 
        WHERE pdb_id LIKE '{pdb_id}' 
        AND chain_type LIKE 'H' 
        AND order_serial >= '{start_pos_para}'
        AND order_serial <= '{end_pos_para}'
        """

        cur.execute(QUERY.format(pdb_id=a[1], start_pos_para=a[4], end_pos_para=a[5]))
        results.append(cur.fetchall())

    file_name = f"energy_distribution_{para_motif_heavy_chain}_{location}.csv"\
        .format(para_motif_heavy_chain=para_motif_heavy_chain, location=location)

    with open(OUT_PATH + file_name, "w") as f:
        out_file = csv.writer(f, lineterminator="\n")
        for r in results:
            out_file.writerow([float(e[0]) for e in r])


    cur.close()
    conn.close()

# get_energy_distribution_for_motif("XXXX", "CDR-H3")
# get_energy_distribution_for_motif("XXXXX", "CDR-H3")
# get_energy_distribution_for_motif("XXXXXX", "CDR-H3")
# get_energy_distribution_for_motif("XXXXXXX", "CDR-H3")
# get_energy_distribution_for_motif("XXXXXXXXX", "CDR-H3")

#get_energy_distribution_for_motif("XXX", "CDR-L1")
#get_energy_distribution_for_motif("XXXX", "LFR3")
#get_energy_distribution_for_motif("XXXX1X", "CDR-H3")
#get_energy_distribution_for_motif("XXXX1XX", "CDR-H3")
#get_energy_distribution_for_motif("XXXXX1X", "CDR-H3")

def get_motifs_pos_for_complex(pdb_id):
    conn = get_conn()
    cur = conn.cursor()

    SQL_QUERY_HC = """
    SELECT para_start_pos, para_end_pos
    FROM motifs
    WHERE pdb_id LIKE '{pdb_id}'
    AND location IN ('CDR-H1', 'CDR-H2', 'CDR-H3', 'HFR1', 'HFR2', 'HFR3', 'HFR4')
    """

    SQL_QUERY_LC = """
    SELECT para_start_pos, para_end_pos
    FROM motifs
    WHERE pdb_id LIKE '{pdb_id}'
    AND location IN ('CDR-L1', 'CDR-L2', 'CDR-L3', 'LFR1', 'LFR2', 'LFR3', 'LFR4')
    """


    cur.execute(SQL_QUERY_HC.format(pdb_id=pdb_id))
    motifs_positions_hc = cur.fetchall()

    cur.execute(SQL_QUERY_LC.format(pdb_id=pdb_id))
    motifs_positions_lc = cur.fetchall()

    cur.close()
    conn.close()

    return motifs_positions_hc, motifs_positions_lc

def get_serial_positions_for_pdb(pdb_id):
    conn = get_conn()
    cur = conn.cursor()

    SQL_QUERY = """
    SELECT MAX(order_serial)
    FROM residues
    WHERE chain_type LIKE '{chain_type}'
    AND pdb_id LIKE '{pdb_id}'
    
    """

    cur.execute(SQL_QUERY.format(chain_type="H", pdb_id=pdb_id))
    max_hc = cur.fetchall()
    cur.execute(SQL_QUERY.format(chain_type="L", pdb_id=pdb_id))
    max_lc = cur.fetchall()

    cur.close()
    conn.close()

    return max_hc, max_lc

def get_ala_ddg_for_pdb(pdb_id):
    from psycopg2.errors import SyntaxError as SyntaxError
    conn = get_conn()
    cur = conn.cursor()

    motifs_positions_hc, motifs_positions_lc = get_motifs_pos_for_complex(pdb_id)
    max_hc, max_lc = get_serial_positions_for_pdb(pdb_id)

    hc_pos_motif = list(itertools.chain(*[list(range(a[0], a[1]+1)) for a in motifs_positions_hc]))
    lc_pos_motif = list(itertools.chain(*[list(range(a[0], a[1] + 1)) for a in motifs_positions_lc]))

    hc_pos_non_motif = [i for i in range(1, max_hc[0][0]+1) if i not in hc_pos_motif]
    lc_pos_non_motif = [i for i in range(1, max_hc[0][0]+1) if i not in lc_pos_motif]

    QUERY = """
    SELECT AVG(ala_scan_energy) 
    FROM residues
    WHERE pdb_id LIKE '{pdb_id}'
    AND chain_type LIKE '{chain_type}'
    AND order_serial IN {pos_list}
    """

    try:
        cur.execute(QUERY.format(pdb_id=pdb_id, chain_type="H", pos_list=tuple(hc_pos_motif)))
        hc_motif_energy = cur.fetchall()
    except SyntaxError:
        conn.rollback()
        hc_motif_energy = [(0.0,)]

    cur.execute(QUERY.format(pdb_id=pdb_id, chain_type="H", pos_list=tuple(hc_pos_non_motif)))
    hc_non_motif_energy = cur.fetchall()

    try:
        cur.execute(QUERY.format(pdb_id=pdb_id, chain_type="L", pos_list=tuple(lc_pos_motif)))
        lc_motif_energy = cur.fetchall()
    except SyntaxError:
        conn.rollback()
        lc_motif_energy = [(0.0,)]

    cur.execute(QUERY.format(pdb_id=pdb_id, chain_type="L", pos_list=tuple(lc_pos_non_motif)))
    lc_non_motif_energy = cur.fetchall()


    cur.close()
    conn.close()

    return hc_motif_energy[0][0], hc_non_motif_energy[0][0], lc_motif_energy[0][0], lc_non_motif_energy[0][0]

def get_avg_energies_for_all():
    conn = get_conn()
    cur = conn.cursor()

    QUERY = """
    SELECT p.pdb_id, i.interaction_energy
    FROM pdbs p 
    INNER JOIN complex_analysis_interaction i
    ON p.pdb_id = i.pdb_id
    """

    cur.execute(QUERY)
    res = cur.fetchall()

    with open(OUT_PATH + "motif_ala_energies.csv", "w") as f:
        out_file = csv.writer(f, lineterminator="\n")
        out_file.writerow(['pdb_id', 'interaction_energy', 'hc_motif_energy', 'hc_non_motif_energy', 'lc_motif_energy', 'lc_non_motif_energy'])
        for r in res:
            pdb_id = r[0]
            print(pdb_id)
            interaction_energy = r[1]
            hc_motif_energy, hc_non_motif_energy, lc_motif_energy, lc_non_motif_energy = get_ala_ddg_for_pdb(pdb_id)
            out_file.writerow([pdb_id, interaction_energy, hc_motif_energy, hc_non_motif_energy, lc_motif_energy, lc_non_motif_energy])



    cur.close()
    conn.close()

get_avg_energies_for_all()