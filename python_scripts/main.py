from Search_and_download_homologues import *
from Process_homologues_and_model import *
from Proteins_score import *

print("--------------------MHAM------------------------\n"
      "------Multi-template Homology Auto Modeling-----\n"
     "-------------------------------------------------\n"
    "-------------------Creators----------------------\n"
    "--------A.Shemyakina--I.Piankov--P.Popov---------\n"
    "-------------------------------------------------\n"
    "This tool is designed to automate the process----\n"
    "of protein modeling.To search for homologues in--\n"
    "-------the Swiss-Model and wwPDB databases-------\n"
    "Alignment is done using mufft with L-INS-i-------\n"
    "algorithm tool.Rosetta is used to model proteins-\n"
    "Ornate 3D is used to assess the quality of models\n"
    "-------------------------------------------------\n"
    "---------------------Note------------------------\n"
    "Protein sequence length is not more than 1000----\n"
    "amino acid residues.-----------------------------\n"
    "Specify only the absolute path to the folder!----\n"
    "-------------------Good luck---------------------\n"
    "-------------------------------------------------\n")

###############################################################################
# Search and download homologues
###############################################################################

print('\n' + '#' * 80 + '\nPart 1: Search and download homologues\n' + '#' * 80 + '\n')

id_protein = input('Target protein UniProt ID (Ex:Q86WT6): ')
directory = input('File save path (Ex:path/to/save/folder/): ')
path_to_Rosetta_dir = input('Provide path to compiled Rosetta directory (Ex:path/to/Rosetta/): ')
path_to_Ornate_dir = input('Provide path to compiled Ornate directory ((Ex:path/to/Ornate/): ')
desired_number_of_models = input('How many models do you want? ')

mean_and_sd_for_function = True
graph_for_function = True

if input("Do you want to calculate the mean and SD for the resulting models? [Y/N]:") in ["n", "N"]:
    mean_and_sd_for_function = False
if input("Do you want to plot with the a quality score for the residuals ? [Y/N]:") in ["n", "N"]:
    graph_for_function = True

if not os.path.exists(directory):
    os.makedirs(directory)

requests_from_uniprot(id_protein, directory)

if input("Do you have a protein database created by this script earlier? [Y/N]:") in ["y", "Y"]:
    path_to_db = input("Path to this database (Ex:path/to/database/database): ")
    print("Start to search homologues using mafft")
    n_in = input('Number homologues (Ex:10):')
    th_in = input('Threshold (Ex:35): ')
    start_mafft(directory, id_protein, n=n_in, th=th_in, flag_path=True, path=path_to_db)
else:
    only_swiss_model = False
    only_wwpdb_model = False
    flag_for_duplicates = True
    if input("Do you want to change the databases used for searching? [Y/N]:") in ["y", "Y"]:
        if input("Use wwPDB database? [Y/N]:") in ["y", "Y"]:
            only_wwpdb_model = True
        elif input("Use Swiss-Model database? [Y/N]:") in ["y", "Y"]:
            only_swiss_model = True
            clean_pdb_db = False
            flag_for_duplicates = False
    if flag_for_duplicates:
        if input("Do you want to clean wwPDB from duplicates? [Y/N]:") in ["y", "Y"]:
            print("Start this process. It may take a few minutes or more.")
            downloading_databases(directory, only_swiss_model, only_wwpdb_model)
        else:
            print("Start this process. It may take a few minutes or more.")
            downloading_databases(directory, only_swiss_model, only_wwpdb_model, clean_pdb_db=False)
    else:
        print("Start this process. It may take a few minutes or more.")
        downloading_databases(directory, only_swiss_model, only_wwpdb_model, clean_pdb_db=False)
    create_database_for_mafft(directory)
    print("Start to search homologues using mafft")
    n = input('Number homologues (Ex:10):')
    th = input('Threshold (Ex:35):')
    start_mafft(directory, id_protein, n, th, flag_path=False)


pdb_download(directory, id_protein)

###############################################################################
# Process homologues and model
###############################################################################

print('\n' + '#' * 80 + '\nPart 2: Process homologues and model\n' + '#' * 80 + '\n')


clean_target_name(directory, id_protein)
make_output_dirs_for_part2(directory)



run_clean_pdb(directory, path_to_Rosetta_dir + 'main/tools/protein_tools/scripts/')
number_of_fastas_in_raw_alignment = mafft_multiple_alignment(directory, id_protein, 'raw_alignment')
number_of_fastas_in_deduplicated_alignment = remove_duplicates('raw_alignment', number_of_fastas_in_raw_alignment, directory, 'raw_alignment_deduplicated')
calculate_identity_for_mult_templates(directory, 'raw_alignment_deduplicated', number_of_fastas_in_deduplicated_alignment)
identity_threshold = input('Provide the identity threshold. Templates with identity below threshold will be removed. ')
if not identity_threshold.isdigit():
    print('Didn\'t catch your answer (not digit), will set identity threshold at 15%')
    identity_threshold = 15
else:
    identity_threshold = int(identity_threshold)
remove_templates_under_threshold(directory, 'identity', identity_threshold)
number_of_fastas_in_clean_alignment = mafft_multiple_alignment(directory, id_protein, 'clean_alignment')
coverage_percent = coverage_minimum(directory, 'clean_alignment', number_of_fastas_in_clean_alignment)
proceed_question = input(f'The coverage of your target protein by selected homologues is {coverage_percent}%. Do you want to proceed? (yes/no): ')
if proceed_question == 'no':
    sys.exit("The execution of the code was interrupted because the homologues do not sufficiently cover the target protein")
elif proceed_question != 'yes':
    print('Didn\'t catch your answer, let\'s say your answer was yes')

convert_fasta_alingment_to_grishin(directory, 'clean_alignment', number_of_fastas_in_clean_alignment)
thread_selected_pdbs(directory, id_protein, path_to_Rosetta_dir)

make_fragments(directory, id_protein, path_to_Rosetta_dir)

modify_xml_file(directory, id_protein)

if not desired_number_of_models.isdigit():
    print('Didn\'t catch your answer (not digit), will make 5 models')
    desired_number_of_models = 5
run_modeling(directory, path_to_Rosetta_dir, id_protein, desired_number_of_models)

###############################################################################
# Scoring
###############################################################################

print('\n' + '#' * 80 + '\nPart 3: Model quality control\n' + '#' * 80 + '\n')


make_output_dirs_for_part3(directory)

run_ornate_score(directory, id_protein, path_to_Ornate_dir)

score_counter(directory, mean_and_sd=mean_and_sd_for_function, graph=graph_for_function)
