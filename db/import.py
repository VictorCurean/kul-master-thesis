import psycopg2
from os import listdir
from os.path import isfile, join

from motifs.structure import Complex
from motifs.utils import AA_3to1

PATH = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\"
path_pdb_files = PATH + "datasets\\"
path_motifs = PATH + "motifs_akbar\\"
path_complex_analysis = PATH + "foldx\\out.complex_analysis\\"
path_ala_scan = PATH + "foldx\\out.ala_scan\\"

DB_NAME = "kul_thesis"
USER = "postgres"
PASSWORD = "admin"
HOST = "localhost"

def get_conn():
    return psycopg2.connect(host=HOST, database=DB_NAME, user=USER, password=PASSWORD)

def test_conn():
    conn = get_conn()
    cur = conn.cursor()
    cur.execute('SELECT version()')
    db_version = cur.fetchone()
    print(db_version)
    cur.close()

def create_tables():
    CREATE_PDB_TABLE = """
    CREATE TABLE IF NOT EXISTS pdbs (
	id serial PRIMARY KEY,
	pdb_id VARCHAR UNIQUE NOT NULL,
	heavy_chain_ann VARCHAR NOT NULL,
	light_chain_ann VARCHAR NOT NULL,
	antigen_ann VARCHAR NOT NULL,
	affinity FLOAT NULL
);
    """

    CREATE_COMPLEX_ANALYSIS_SUMMARY_TABLE = """
    CREATE TABLE IF NOT EXISTS complex_analysis_summary (
    id serial PRIMARY KEY,
    pdb_id VARCHAR UNIQUE NOT NULL,
    intraclashes_ab FLOAT NOT NULL,
    intraclashes_ag FLOAT NOT NULL,
    total_energy FLOAT NOT NULL,
    stability_ab FLOAT NOT NULL,
    stability_ag FLOAT NOT NULL
    );
    """

    CREATE_COMPLEX_ANALYSIS_INDIVIDUAL_ENERGIES_TABLE = """
    CREATE TABLE IF NOT EXISTS complex_analysis_individual_energies (
    id serial PRIMARY KEY,
    pdb_id VARCHAR UNIQUE NOT NULL,
    
    total_energy_ab FLOAT NOT NULL,
    backbone_hbond_ab FLOAT NOT NULL,
    sidechain_hbond_ab FLOAT NOT NULL,
    vdw_ab FLOAT NOT NULL,
    electrostatics_ab FLOAT NOT NULL,
    solvation_polar_ab FLOAT NOT NULL,
    solvation_hydrophobic_ab FLOAT NOT NULL,
    vdw_clashes_ab FLOAT NOT NULL,
    entropy_sidechain_ab FLOAT NOT NULL,
    entropy_mainchain_ab FLOAT NOT NULL,
    sloop_entropy_ab FLOAT NOT NULL,
    mloop_entropy_ab FLOAT NOT NULL,
    cis_bond_ab FLOAT NOT NULL,
    torsional_clash_ab FLOAT NOT NULL,
    backbone_clash_ab FLOAT NOT NULL,
    helix_dipole_ab FLOAT NOT NULL,
    water_bridge_ab FLOAT NOT NULL,
    disulfide_ab FLOAT NOT NULL,
    electrostatic_kon_ab FLOAT NOT NULL,
    partial_covalent_bonds_ab FLOAT NOT NULL,
    energy_ionisation_ab FLOAT NOT NULL,
    entropy_complex_ab FLOAT NOT NULL,
    
    total_energy_ag FLOAT NOT NULL,
    backbone_hbond_ag FLOAT NOT NULL,
    sidechain_hbond_ag FLOAT NOT NULL,
    vdw_ag FLOAT NOT NULL,
    electrostatics_ag FLOAT NOT NULL,
    solvation_polar_ag FLOAT NOT NULL,
    solvation_hydrophobic_ag FLOAT NOT NULL,
    vdw_clashes_ag FLOAT NOT NULL,
    entropy_sidechain_ag FLOAT NOT NULL,
    entropy_mainchain_ag FLOAT NOT NULL,
    sloop_entropy_ag FLOAT NOT NULL,
    mloop_entropy_ag FLOAT NOT NULL,
    cis_bond_ag FLOAT NOT NULL,
    torsional_clash_ag FLOAT NOT NULL,
    backbone_clash_ag FLOAT NOT NULL,
    helix_dipole_ag FLOAT NOT NULL,
    water_bridge_ag FLOAT NOT NULL,
    disulfide_ag FLOAT NOT NULL,
    electrostatic_kon_ag FLOAT NOT NULL,
    partial_covalent_bonds_ag FLOAT NOT NULL,
    energy_ionisation_ag FLOAT NOT NULL,
    entropy_complex_ag FLOAT NOT NULL
    );
    """

    CREATE_COMPLEX_ANALYSIS_INTERACTION_TABLE = """
    CREATE TABLE IF NOT EXISTS complex_analysis_interaction (
    id SERIAL PRIMARY KEY,
    pdb_id VARCHAR UNIQUE NOT NULL,
    
    interaction_energy FLOAT NOT NULL,
    backbone_hbond FLOAT NOT NULL,
    sidechain_hbond FLOAT NOT NULL,
    vdw FLOAT NOT NULL,
    electrostatics FLOAT NOT NULL,
    solvation_polar FLOAT NOT NULL,
    solvation_hydrophobic FLOAT NOT NULL,
    vdw_clashes FLOAT NOT NULL,
    entropy_sidechain FLOAT NOT NULL, 
    entropy_mainchain FLOAT NOT NULL,
    sloop_entropy FLOAT NOT NULL,
    mloop_entropy FLOAT NOT NULL,
    cis_bond FLOAT NOT NULL,
    torsional_clash FLOAT NOT NULL,
    backbone_clash FLOAT NOT NULL,
    helix_dipole FLOAT NOT NULL,
    water_bridge FLOAT NOT NULL,
    disulfide FLOAT NOT NULL,
    electrostatic_kon FLOAT NOT NULL,
    partial_covalent_bonds FLOAT NOT NULL,
    energy_ionisation FLOAT NOT NULL,
    entropy_complex FLOAT NOT NULL,
    number_of_residues INTEGER NOT NULL,
    interface_residues INTEGER NOT NULL,
    interface_residues_clashing INTEGER NOT NULL,
    interface_residues_vdw_clashing INTEGER NOT NULL,
    interface_residues_bb_clashing INTEGER NOT NULL
    )
    """

    CREATE_RESIDUE_TABLE = """
    CREATE TABLE IF NOT EXISTS residues (
    id SERIAL PRIMARY KEY,
    pdb_id VARCHAR NOT NULL,
    chain_type VARCHAR NOT NULL,
    amino_acid VARCHAR NOT NULL,
    order_serial INTEGER NOT NULL,
    order_martin INTEGER NOT NULL,
    ala_scan_energy FLOAT NOT NULL
    )
    """

    CREATE_MOTIFS_TABLE = """
    CREATE TABLE IF NOT EXISTS motifs (
    id SERIAL PRIMARY KEY,
    pdb_id VARCHAR NOT NULL,
    para_motif VARCHAR NOT NULL,
    epi_motif VARCHAR NOT NULL,
    para_start_pos INTEGER NOT NULL,
    para_end_pos INTEGER NOT NULL,
    epi_start_pos INTEGER NOT NULL,
    epi_end_pos INTEGER NOT NULL,
    interacting_residues_para INTEGER NOT NULL,
    interacting_residues_epi INTEGER NOT NULL,
    location VARCHAR NOT NULL
    )
    """

    conn = get_conn()
    cur = conn.cursor()

    # cur.execute(CREATE_PDB_TABLE)
    # cur.execute(CREATE_COMPLEX_ANALYSIS_SUMMARY_TABLE)
    # cur.execute(CREATE_COMPLEX_ANALYSIS_INDIVIDUAL_ENERGIES_TABLE)
    # cur.execute(CREATE_COMPLEX_ANALYSIS_INTERACTION_TABLE)
    # cur.execute(CREATE_RESIDUE_TABLE)
    cur.execute(CREATE_MOTIFS_TABLE)

    conn.commit()
    cur.close()

def get_ag_ann(file):
    with open(file, "r") as file:
        for i, line in enumerate(file):
            if i == 8:
                return line.split("    ")[0].split(" ")[3]

def insert_complex_analysis_summary(pdb_id, conn):
    file_summary = path_complex_analysis + "Summary_" + pdb_id + "_AC.fxout"
    vals = list()

    with open(file_summary, "r") as f:
        for i, line in enumerate(f):
            if i == 9:
                vals = line.split("\t")
                # intraclashes_ab = line.split("\t")[3]
                # intraclashes_ag = line.split("\t")[4]
                # total_energy = line.split("\t")[5]
                # stability_ab = line.split("\t")[6]
                # stability_ag = line.split("\t")[7]

    cur = conn.cursor()
    INSERT_SQL = """
    INSERT INTO complex_analysis_summary(pdb_id, intraclashes_ab, intraclashes_ag, total_energy, stability_ab, stability_ag)
    VALUES('{pdb_id}', '{intraclashes_ab}', '{intraclashes_ag}', '{total_energy}', '{stability_ab}', '{stability_ag}')
    """
    cur.execute(INSERT_SQL.format(pdb_id=pdb_id, intraclashes_ab=vals[3], intraclashes_ag=vals[4],
                                  total_energy=vals[5], stability_ab=vals[6], stability_ag=vals[7]))

    conn.commit()
    cur.close()

def insert_complex_analysis_indiv_energies(pdb_id, conn):
    file_indiv_energies = path_complex_analysis + "Indiv_energies_" + pdb_id + "_AC.fxout"
    vals_ab = list()
    vals_ag = list()

    cur = conn.cursor()

    with open(file_indiv_energies, "r") as f:
        for i, line in enumerate(f):
            if i == 9:
                #Antibody
                vals_ab = line.split("\t")

            elif i == 10:
                #Antigen
                vals_ag = line.split("\t")


    INSERT_SQL = """
    INSERT INTO complex_analysis_individual_energies(pdb_id, total_energy_ab, backbone_hbond_ab, sidechain_hbond_ab, vdw_ab, electrostatics_ab, solvation_polar_ab, solvation_hydrophobic_ab, vdw_clashes_ab, entropy_sidechain_ab, entropy_mainchain_ab, sloop_entropy_ab, mloop_entropy_ab, cis_bond_ab, torsional_clash_ab, backbone_clash_ab, helix_dipole_ab, water_bridge_ab, disulfide_ab, electrostatic_kon_ab, partial_covalent_bonds_ab, energy_ionisation_ab, entropy_complex_ab, 
    total_energy_ag, backbone_hbond_ag, sidechain_hbond_ag, vdw_ag, electrostatics_ag, solvation_polar_ag, solvation_hydrophobic_ag, vdw_clashes_ag, entropy_sidechain_ag, entropy_mainchain_ag, sloop_entropy_ag, mloop_entropy_ag, cis_bond_ag, torsional_clash_ag, backbone_clash_ag, helix_dipole_ag, water_bridge_ag, disulfide_ag, electrostatic_kon_ag, partial_covalent_bonds_ag, energy_ionisation_ag, entropy_complex_ag)
    VALUES ('{pdb_id}', '{total_energy_ab}', '{backbone_hbond_ab}', '{sidechain_hbond_ab}', '{vdw_ab}', '{electrostatics_ab}', '{solvation_polar_ab}', '{solvation_hydrophobic_ab}', '{vdw_clashes_ab}', '{entropy_sidechain_ab}', '{entropy_mainchain_ab}', '{sloop_entropy_ab}', '{mloop_entropy_ab}', '{cis_bond_ab}', '{torsional_clash_ab}', '{backbone_clash_ab}', '{helix_dipole_ab}', '{water_bridge_ab}', '{disulfide_ab}', '{electrostatic_kon_ab}', '{partial_covalent_bonds_ab}', '{energy_ionisation_ab}', '{entropy_complex_ab}', 
    '{total_energy_ag}', '{backbone_hbond_ag}', '{sidechain_hbond_ag}', '{vdw_ag}', '{electrostatics_ag}', '{solvation_polar_ag}', '{solvation_hydrophobic_ag}', '{vdw_clashes_ag}', '{entropy_sidechain_ag}', '{entropy_mainchain_ag}', '{sloop_entropy_ag}', '{mloop_entropy_ag}', '{cis_bond_ag}', '{torsional_clash_ag}', '{backbone_clash_ag}', '{helix_dipole_ag}', '{water_bridge_ag}', '{disulfide_ag}', '{electrostatic_kon_ag}', '{partial_covalent_bonds_ag}', '{energy_ionisation_ag}', '{entropy_complex_ag}'
    )
    """

    cur.execute(INSERT_SQL.format(pdb_id=pdb_id, total_energy_ab=vals_ab[2], backbone_hbond_ab=vals_ab[3], sidechain_hbond_ab=vals_ab[4], vdw_ab=vals_ab[5], electrostatics_ab=vals_ab[6], solvation_polar_ab=vals_ab[7], solvation_hydrophobic_ab=vals_ab[8], vdw_clashes_ab=vals_ab[9], entropy_sidechain_ab=vals_ab[10], entropy_mainchain_ab=vals_ab[11], sloop_entropy_ab=vals_ab[12], mloop_entropy_ab=vals_ab[13], cis_bond_ab=vals_ab[14], torsional_clash_ab=vals_ab[15], backbone_clash_ab=vals_ab[16], helix_dipole_ab=vals_ab[17], water_bridge_ab=vals_ab[18], disulfide_ab=vals_ab[19], electrostatic_kon_ab=vals_ab[20], partial_covalent_bonds_ab=vals_ab[21], energy_ionisation_ab=vals_ab[22], entropy_complex_ab=vals_ab[23],
                                  total_energy_ag=vals_ag[2], backbone_hbond_ag=vals_ag[3], sidechain_hbond_ag=vals_ag[4], vdw_ag=vals_ag[5], electrostatics_ag=vals_ag[6], solvation_polar_ag=vals_ag[7], solvation_hydrophobic_ag=vals_ag[8], vdw_clashes_ag=vals_ag[9], entropy_sidechain_ag=vals_ag[10], entropy_mainchain_ag=vals_ag[11], sloop_entropy_ag=vals_ag[12], mloop_entropy_ag=vals_ag[13], cis_bond_ag=vals_ag[14], torsional_clash_ag=vals_ag[15], backbone_clash_ag=vals_ag[16], helix_dipole_ag=vals_ag[17], water_bridge_ag=vals_ag[18], disulfide_ag=vals_ag[19], electrostatic_kon_ag=vals_ag[20], partial_covalent_bonds_ag=vals_ag[21], energy_ionisation_ag=vals_ag[22], entropy_complex_ag=vals_ag[23]))

    conn.commit()
    cur.close()

def insert_complex_analysis_interaction(pdb_id, conn):
    file_interaction_energies = path_complex_analysis + "Interaction_" + pdb_id + "_AC.fxout"
    vals = list()

    cur = conn.cursor()

    with open(file_interaction_energies, "r") as f:
        for i, line in enumerate(f):
            if i == 9:
                vals = line.split("\t")

    INSERT_SQL = """
    INSERT INTO complex_analysis_interaction(pdb_id, interaction_energy, backbone_hbond, sidechain_hbond, vdw, electrostatics, solvation_polar, solvation_hydrophobic, vdw_clashes, entropy_sidechain, entropy_mainchain, sloop_entropy, mloop_entropy, cis_bond, torsional_clash, backbone_clash, helix_dipole, water_bridge, disulfide, electrostatic_kon, partial_covalent_bonds, energy_ionisation, entropy_complex, number_of_residues, interface_residues, interface_residues_clashing, interface_residues_vdw_clashing, interface_residues_bb_clashing)
    VALUES ('{pdb_id}', '{interaction_energy}', '{backbone_hbond}', '{sidechain_hbond}', '{vdw}', '{electrostatics}', '{solvation_polar}', '{solvation_hydrophobic}', '{vdw_clashes}', '{entropy_sidechain}', '{entropy_mainchain}', '{sloop_entropy}', '{mloop_entropy}', '{cis_bond}', '{torsional_clash}', '{backbone_clash}', '{helix_dipole}', '{water_bridge}', '{disulfide}', '{electrostatic_kon}', '{partial_covalent_bonds}', '{energy_ionisation}', '{entropy_complex}', '{number_of_residues}', '{interface_residues}', '{interface_residues_clashing}', '{interface_residues_vdw_clashing}', '{interface_residues_bb_clashing}')
    """

    cur.execute(INSERT_SQL.format(pdb_id=pdb_id, interaction_energy=vals[5], backbone_hbond=vals[6], sidechain_hbond=vals[7], vdw=vals[8], electrostatics=vals[9], solvation_polar=vals[10], solvation_hydrophobic=vals[11], vdw_clashes=vals[12], entropy_sidechain=vals[13], entropy_mainchain=vals[14], sloop_entropy=vals[15], mloop_entropy=vals[16], cis_bond=vals[17], torsional_clash=vals[18], backbone_clash=vals[19], helix_dipole=vals[20], water_bridge=vals[21], disulfide=vals[22], electrostatic_kon=vals[23], partial_covalent_bonds=vals[24], energy_ionisation=vals[25], entropy_complex=vals[26], number_of_residues=vals[27], interface_residues=vals[28], interface_residues_clashing=vals[29], interface_residues_vdw_clashing=vals[30], interface_residues_bb_clashing=vals[31]))

    conn.commit()
    cur.close()

def insert_residue_information(pdb_id, conn):

    cur = conn.cursor()

    #order - L, H, AG
    ala_scan_file = path_ala_scan + pdb_id + "_AS.fxout"
    order = ["L", "H", "AG", "AG", "AG"]
    order_index = 0


    with open(ala_scan_file, "r") as f:
        nr_prev = 0
        serial_numbering = 1

        for i, line in enumerate(f):
            vals = line.split(" ")
            aa_martin = int(vals[1])

            if aa_martin < nr_prev:
                order_index += 1
                serial_numbering = 1

            aa_3_not = vals[0]
            energy_change = float(vals[7].replace("\n", ""))
            aa_1_not = AA_3to1[aa_3_not]
            chain_type = order[order_index]

            INSERT_SQL = """
            INSERT INTO residues(pdb_id, chain_type, amino_acid, order_serial, order_martin, ala_scan_energy)
            VALUES ('{pdb_id}', '{chain_type}', '{amino_acid}', '{order_serial}', '{order_martin}', '{ala_scan_energy}')
            """

            cur.execute(INSERT_SQL.format(pdb_id=pdb_id, chain_type=chain_type, amino_acid=aa_1_not, order_serial=serial_numbering, order_martin=aa_martin, ala_scan_energy=energy_change))

            nr_prev = aa_martin
            serial_numbering += 1

    conn.commit()
    cur.close()

def insert_motifs():
    path = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\datasets"
    files = [f for f in listdir(path) if isfile(join(path, f))]

    conn = get_conn()
    cur = conn.cursor()

    for f in files:
        print(f)
        comp = Complex(f, path + "\\" + f)
        comp.get_motifs()

        pdb_id = f.split(".")[0]

        for motif in comp.motifs:
            para_motif = motif.get_akbar_notation_from_res(motif.encodings_origin)
            epi_motif = motif.get_akbar_notation_from_res(motif.encodings_target)
            interacting_residues_para = para_motif.count("X")
            interacting_residues_epi = epi_motif.count("X")
            start_pos_para = min(motif.positions_ab)
            end_pos_para = max(motif.positions_ab)
            start_pos_epi = min(motif.positions_ag)
            end_pos_epi = max(motif.positions_ag)
            motif_para = motif.get_akbar_notation_from_res(motif.encodings_origin)
            motif_epi = motif.get_akbar_notation_from_res(motif.encodings_target)
            location = motif.ab_location

            INSERT_SQL = """INSERT INTO motifs(pdb_id, para_motif, epi_motif, para_start_pos, para_end_pos, epi_start_pos, epi_end_pos, interacting_residues_para, interacting_residues_epi, location)
            VALUES ('{pdb_id}', '{para_motif}', '{epi_motif}', '{para_start_pos}', '{para_end_pos}', '{epi_start_pos}', '{epi_end_pos}', '{interacting_residues_para}', '{interacting_residues_epi}', '{location}')"""

            cur.execute(INSERT_SQL.format(pdb_id=pdb_id, para_motif=motif_para, epi_motif=motif_epi, para_start_pos=start_pos_para, para_end_pos=end_pos_para, epi_end_pos=end_pos_epi, epi_start_pos=start_pos_epi, interacting_residues_para=interacting_residues_para, interacting_residues_epi=interacting_residues_epi, location=location))
            conn.commit()
    cur.close()
    conn.close()






def insert_pdbs():
    path = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\datasets\\"
    files = [f for f in listdir(path) if isfile(join(path, f))]

    conn = get_conn()
    cur = conn.cursor()

    for f in files:
        pdb_id = f.split(".")[0]
        hc_ann = "H"
        lc_ann = "L"
        ag_ann = get_ag_ann(path + f)

        INSERT_SQL = """
        INSERT INTO pdbs(pdb_id, heavy_chain_ann, light_chain_ann, antigen_ann)
        VALUES ('{pdb_id}', '{heavy_chain_ann}', '{light_chain_ann}', '{antigen_ann}')
        """
        cur.execute(INSERT_SQL.format(pdb_id=pdb_id, heavy_chain_ann=hc_ann, light_chain_ann=lc_ann, antigen_ann=ag_ann))

        insert_complex_analysis_summary(pdb_id, conn)
        insert_complex_analysis_indiv_energies(pdb_id, conn)
        insert_complex_analysis_interaction(pdb_id, conn)
        insert_residue_information(pdb_id, conn)

        conn.commit()

    cur.close()
    conn.close()

create_tables()
#insert_pdbs()

insert_motifs()