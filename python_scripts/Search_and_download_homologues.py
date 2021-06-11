import requests
from Bio import SeqIO
import subprocess
import urllib.request
import sys
import os
import os.path


def requests_from_uniprot(uniprot_id, directory):
    '''
    Function takes Uniprot ID \n
    Return file.fasta of this target protein

    :param uniprot_id: id of target protein
    :type uniprot_id: str
    :param directory: dir  where the fasta file  will be saved
    :type directory: str
    :return: fasta file with target protein
    :rtype: fasta.file
    '''
    base = 'http://www.uniprot.org'
    kb_endpoint = '/uniprot/'
    payload = {'query': f'{uniprot_id}',
               'format': 'fasta'}
    result = requests.get(base + kb_endpoint, params=payload)
    if result.ok:
        with open(os.path.join(directory, f"{uniprot_id}.fasta"), "w") as f:
            f.write(result.text)
    else:
        return print('Something went wrong: ', result.status_code)


def clean_db_pdb(name, directory):
    '''
    This function cleaning wwPDB database from duplicates

    :param name: name fasta file
    :type name: str
    :param directory: path to folder where will be saved file
    :type directory: str
    :return: clear database
    :rtype: file
    '''
    fasta_sequences = SeqIO.parse(open(name), 'fasta')
    first_elem = str()
    with open(os.path.join(directory, "clean_pdb_db.fasta"), "w") as f:
        for fasta in fasta_sequences:
            seq = fasta.seq
            if seq != first_elem:
                first_elem = seq
                SeqIO.write(fasta, f, 'fasta')


def downloading_databases(directory, only_swiss_model=False, only_wwpdb_model=False, clean_pdb_db=True):
    '''
    Function take path to directory and flag for download db. In this function user can choose which db will be used
    for search homologous. The flags are used to select the db.

    :param directory: dir where the db file will be saved
    :type directory: str
    :param only_swiss_model: if choose this flag only swiss_model.db will be download and will be used for search of
    homologous
    :type only_swiss_model: bool
    :param only_wwpdb_model: if choose this flag wwpdb.db will be download and will be used for search of
    homologous
    :type only_wwpdb_model: bool
    :param clean_pdb_db: if choose this flag wwpdb.db will be cleaned from duplicates with using function clean_db_pdb
    :type clean_pdb_db: bool
    :return: file with merge db or not merge db
    :rtype: file
    '''
    out_uniprot = os.path.join(directory, "uniprot_sprot.fasta.gz")
    out_pdb_db = os.path.join(directory, "pdb_seqres.fasta")
    new_name_out_uniprot = os.path.join(directory, "uniprot_sprot.fasta")
    out_protein_db = os.path.join(directory, "db_for_protein.fasta")
    if only_swiss_model:
        urllib.request.urlretrieve(
            "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
            out_uniprot)
        gunzip_uniprot = subprocess.Popen([f"gunzip {out_uniprot}"], shell=True, encoding='utf-8')
        gunzip_uniprot = gunzip_uniprot.wait()
        subprocess.Popen([f"mv {new_name_out_uniprot} {out_protein_db}"], shell=True, encoding='utf-8')
    elif only_wwpdb_model:
        urllib.request.urlretrieve("https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt", out_pdb_db)
        if clean_pdb_db:
            clean_db_pdb(out_pdb_db, directory)
            new_name_pdb_db = os.path.join(directory, "clean_pdb_db.fasta")
            subprocess.Popen([f"mv {new_name_pdb_db} {out_protein_db}"], shell=True, encoding='utf-8')
        else:
            subprocess.Popen([f"mv {out_protein_db} {out_protein_db}"], shell=True, encoding='utf-8')
    else:
        urllib.request.urlretrieve(
            "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
            out_uniprot)
        gunzip_uniprot = subprocess.Popen([f"gunzip {out_uniprot}"], shell=True, encoding='utf-8')
        gunzip_uniprot = gunzip_uniprot.wait()
        urllib.request.urlretrieve("https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt", out_pdb_db)
        if clean_pdb_db:
            clean_db_pdb(out_pdb_db, directory)
            new_name_db = os.path.join(directory, "clean_pdb_db.fasta")
            out_protein_db = os.path.join(directory, "db_for_protein.fasta")
            subprocess.Popen([f"cat {new_name_db} {new_name_out_uniprot} > {out_protein_db}"],
                             shell=True, encoding='utf-8')
        else:
            subprocess.Popen([f"cat {out_pdb_db} {new_name_out_uniprot} > {out_protein_db}"],
                             shell=True, encoding='utf-8')
    print("Databases have been downloaded !")


def create_database_for_mafft(directory):
    '''
    Created special database for mafft tools

    :param directory: path to folder
    :type directory: str
    :return: folder with database for mafft
    :rtype: dir
    '''
    name_of_dir_with_db = os.path.join(directory, "db_for_mafft")
    name_protein_db = os.path.join(directory, "db_for_protein.fasta")
    if not os.path.exists(name_of_dir_with_db):
        os.makedirs(name_of_dir_with_db)
    mafft_database = subprocess.Popen(
        [
            f"makeblastdb -in {name_protein_db} -out {name_of_dir_with_db}/db_for_mafft -dbtype prot -parse_seqids -hash_index"],
        shell=True, encoding='utf-8')
    (output, err) = mafft_database.communicate()
    mafft_database = mafft_database.wait()
    if err:
        print(err)
    print("Database for mafft created")


def start_mafft(directory, id_protein, n=30, th=40, flag_path=False, path=None):
    '''
    Start mafft with L-INS-i algorithm
    (more information here:https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html)
    User can choose number of homologoues and level of threshold line.

    :param directory: dir where the db file will be saved
    :type directory: str
    :param id_protein: fasta file with target protein
    :type id_protein: file
    :param n: number of homologoues
    :type n: int
    :param th: of threshold line
    :type th: int
    :param flag_path: if user have db which was created earlier takes this flag
    :type flag_path: bool
    :param path: path to db which was created earlier
    :type path: bool
    :return: file_output.fasta with list of aligned homologues
    :rtype: file
    '''
    name_of_dir_with_db = os.path.join(directory, "db_for_mafft/db_for_mafft")
    path_to_fasta = os.path.join(directory, f"{id_protein}.fasta")
    path_to_out_fasta = os.path.join(directory, f"{id_protein}_output.fasta")
    if flag_path == False:
        path = name_of_dir_with_db
    mafft_homologs = subprocess.Popen([f"mafft-homologs -l -d {path} \
                                        -o '--thread 8 --threadtb 5 --threadit 0 --reorder\
                                        --anysymbol --maxiterate 1000 --retree 1 --localpair' \
                                        -a {n} -e 1.0e-{th} -f {path_to_fasta} > {path_to_out_fasta}"],
                                      shell=True,
                                      encoding='utf-8')
    (output, err) = mafft_homologs.communicate()
    if err:
        print(err)
    mafft_homologs = mafft_homologs.wait()
    print(f"Mafft searched {n} homologues")


def fasta_parser(directory, name_fasta, name_with_chain=False):
    '''
    Easy label of protein fasta parser

    :param directory: dir where the db file will be saved
    :type directory: str
    :param name_fasta: file_output.fasta with list of aligned homologues
    :type name_fasta: file
    :param name_with_chain: choose this flag if user want to save list with name which not only name but name with chain
    :type name_with_chain: bool
    :return: list with id name of protein homologues
    :rtype: list
    '''
    fasta_name_and_path = os.path.join(directory, f"{name_fasta}_output.fasta")
    fasta_sequences = SeqIO.parse(open(fasta_name_and_path), 'fasta')
    for fasta in fasta_sequences:
        break
    name_list = []
    for fasta in fasta_sequences:
        if name_with_chain:
            fasta_name = fasta.id.split("|")
            if fasta_name[1][0].isdigit():
                fasta_name = fasta_name[1] + "_" + fasta_name[2]
            else:
                fasta_name = fasta_name[1]
            name_list.append(fasta_name)
        else:
            name_list.append(fasta.id.split("|")[1])
    return name_list


def pdb_download(directory, name_fasta):
    '''
    Download all homologues in pdb format \n
    Returns files.pdb

    :param directory: dir where the db file will be saved
    :type directory: str
    :param name_fasta: file_output.fasta with list of aligned homologues
    :type name_fasta: file
    :return: dir with downloaded pdb structure
    :rtype: file
    '''
    path_to_orig_templates = os.path.join(directory, "orig_templates")
    if not os.path.exists(path_to_orig_templates):
        os.makedirs(path_to_orig_templates)
    list_pdb_structure = fasta_parser(directory, name_fasta)
    list_with_name = fasta_parser(directory, name_fasta, name_with_chain=True)
    for name in range(0,len(list_pdb_structure)):
        if list_pdb_structure[name][0].isdigit():
            pdbfn = list_pdb_structure[name] + ".pdb"
            name_with_chain = list_with_name[name] + ".pdb"
            url = "https://files.rcsb.org/download/" + pdbfn
            outfnm = os.path.join(path_to_orig_templates, name_with_chain)
            try:
                urllib.request.urlretrieve(url, outfnm)
                print(f"Structure {name_with_chain} successfully downloaded")
            except Exception as err:
                print(str(err), file=sys.stderr)
        else:
            uniname = list_pdb_structure[name] + ".pdb"
            url_uni = "https://swissmodel.expasy.org/repository/uniprot/" + uniname
            outfnm_uni = os.path.join(path_to_orig_templates, uniname)
            try:
                urllib.request.urlretrieve(url_uni, outfnm_uni)
                print(f"Structure {uniname} successfully downloaded")
            except Exception as err:
                print(str(err), file=sys.stderr)
