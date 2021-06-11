from Preparation import *
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
# Preparation
###############################################################################



# Checking dependencies - Ornate and Rosetta

path_to_our_script = os.getcwd() + '/'
if os.path.exists(path_to_our_script + 'paths'):
    with open(path_to_our_script + 'paths') as paths:
        path_to_Ornate_dir = paths.readline().strip()
        print(f'Using Ornate from {path_to_Ornate_dir}')
        path_to_Rosetta_dir = paths.readline().strip()
        print(f'Using Rosetta from {path_to_Rosetta_dir}')
else:
    print('Checking dependencies...')
    path_to_Ornate_dir = find_path_to_dependency(path_to_our_script, 'Ornate', 'score.py')
    path_to_Rosetta_dir = find_path_to_dependency(path_to_our_script, 'Rosetta', 'main')

    if path_to_Ornate_dir is None:
        sys.exit(
            "The script requires Ornate software which is missing on your computer "
            "(we were searching for Ornate directory with score.py file inside it). "
            "Check out https://team.inria.fr/nano-d/software/Ornate/")
    else:
        print(f'Using Ornate from {path_to_Ornate_dir}')

    if path_to_Rosetta_dir is None:
        sys.exit(
            "The script requires Rosetta software which is missing on your computer"
            "(we were searching for Rosetta directory with main directory inside it). "
            "Check out https://www.rosettacommons.org/demos/latest/tutorials/install_build/install_build")
    else:
        print(f'Using Rosetta from {path_to_Rosetta_dir}')

# Asking parameters for script execution from user

print('\nSetting script general parameters')
id_protein, directory = general_questions()
print('\nSetting parameters for protein databases downloading')
path_to_db, only_swiss_model, only_wwpdb_model, clean_pdb_db = questions_about_protein_databases()
print('\nSetting parameters for homologues search')
n_of_homologs, e_value_threshold = questions_about_homologs_search()
print('\nSetting parameters for homologue processing and modelling')
identity_threshold, coverage_threshold, desired_number_of_models = questions_about_homologs_processing_and_modeling()
print('\nSetting parameters for scoring')
mean_and_sd_for_function, graph_for_function = questions_about_scoring()

###############################################################################
# Search and download homologues
###############################################################################

print('\n' + '#' * 80 + '\nPart 1: Search and download homologues\n' + '#' * 80 + '\n')

# Initial preparations for part 1

if not os.path.exists(directory):
    os.makedirs(directory)

# download target sequence

requests_from_uniprot(id_protein, directory)

# download database

if path_to_db:
    print("Start searching homologues using mafft")
    start_mafft(directory, id_protein, n=n_of_homologs, th=e_value_threshold, flag_path=True, path=path_to_db)
else:
    print("Start downloading database(s). It may take a few minutes or more.")
    downloading_databases(directory, only_swiss_model, only_wwpdb_model, clean_pdb_db=clean_pdb_db)
    create_database_for_mafft(directory)
    print("Start searching homologues using mafft")
    start_mafft(directory, id_protein, n=n_of_homologs, th=e_value_threshold, flag_path=False)

# download selected homologues

pdb_download(directory, id_protein)
if len(next(os.walk(directory + 'orig_templates/'))[2]) == 0:
    sys.exit("No homologues were found")

###############################################################################
# Process homologues and model
###############################################################################

print('\n' + '#' * 80 + '\nPart 2: Process homologues and model\n' + '#' * 80 + '\n')

# Initial preparations for part 2

clean_target_name(directory, id_protein)
make_output_dirs_for_part2(directory)

# Processing of homologs

run_clean_pdb(directory, path_to_Rosetta_dir + 'main/tools/protein_tools/scripts/')
number_of_fastas_in_raw_alignment = mafft_multiple_alignment(directory, id_protein, 'raw_alignment')
number_of_fastas_in_deduplicated_alignment = remove_duplicates('raw_alignment', number_of_fastas_in_raw_alignment,
                                                               directory, 'raw_alignment_deduplicated')

# Check identity and coverage percents

calculate_identity_for_mult_templates(directory, 'raw_alignment_deduplicated',
                                      number_of_fastas_in_deduplicated_alignment)
remove_templates_under_threshold(directory, 'identity', identity_threshold)
if len(next(os.walk(directory + 'Modeling/cleaned_template_pdbs/'))[2]) == 0:
    sys.exit(
        "No homologue has passed the threshold of identity")
number_of_fastas_in_clean_alignment = mafft_multiple_alignment(directory, id_protein, 'clean_alignment')
coverage_percent = calculate_coverage(directory, 'clean_alignment', number_of_fastas_in_clean_alignment)
if coverage_percent < coverage_threshold:
    sys.exit(
        f"Homologues do not sufficiently cover the target protein. The coverage percent is {coverage_percent}")

# Modelling

convert_fasta_alingment_to_grishin(directory, 'clean_alignment', number_of_fastas_in_clean_alignment)
thread_selected_pdbs(directory, id_protein, path_to_Rosetta_dir)
# make_fragments(directory, id_protein, path_to_Rosetta_dir)
modify_xml_file(directory, id_protein)
run_modeling(directory, path_to_Rosetta_dir, id_protein, desired_number_of_models)

###############################################################################
# Scoring
###############################################################################

print('\n' + '#' * 80 + '\nPart 3: Model quality control\n' + '#' * 80 + '\n')

make_output_dirs_for_part3(directory)
run_ornate_score(directory, id_protein, path_to_Ornate_dir)
score_counter(directory, mean_and_sd=mean_and_sd_for_function, graph=graph_for_function)
