import os
import subprocess
import shutil
import re
from lxml import etree


def parse_multifasta_file(file, number_of_fastas):
    """
    Analogue of biopython SeqIO.parse (allows to avoid using excess dependencies)

    :param file: path to file containing multiple fasta sequences (ex: path/to/file.fasta)
    :type file: str
    :param number_of_fastas: number of fasta sequences in file
    :type number_of_fastas: int
    :return: returns name (comment line) and sequence of fasta when called with next()
    :rtype: Generator[str, str]
    """

    with open(file) as file:
        for i in range(number_of_fastas):
            fasts_seq = ''
            fasta_name = file.readline().strip()[1:]
            end_of_file = False
            end_of_seq = False
            while not end_of_seq and not end_of_file:
                x = file.tell()
                seq = file.readline()
                if not seq:
                    end_of_file = True
                elif '>' not in seq:
                    fasts_seq = fasts_seq + seq
                else:
                    file.seek(x)
                    end_of_seq = True
            fasts_seq = re.sub(r'\n', '', fasts_seq)
            yield fasta_name, fasts_seq


def clean_target_name(path, target):
    """
    Script requires that in fasta file of target protein comment line contains only its Uniprot ID
    This function overwrites the fasta file simplifying comment line (the file name is preserved) \n
    Returns nothing

    :param path: path to directory with the fasta file (ex: path/to/directory/)
    :type path: str
    :param target: Uniprot ID of target protein (ex: Q86WT6)
    :type target: str
    """

    path_to_file = path + target + '.fasta'
    fasta = parse_multifasta_file(path_to_file, 1)
    target_name, target_seq = next(fasta)
    with open(path_to_file, 'w') as clean_file:
        clean_file.write(f'>{target}\n')
        clean_file.write(target_seq + '\n')


def make_output_dirs_for_part2(parent_dir):
    """
    Functions creates all directories necessary for the second part 'Process homologues and model'
    in provided working directory (takes into account the case if the directory exists already) \n
    Returns nothing

    :param parent_dir: path to working directory (ex: path/to/directory/)
    :type parent_dir: str
    """

    if not os.path.exists(parent_dir + 'Modeling/'):
        os.makedirs(parent_dir + 'Modeling/')
    if not os.path.exists(parent_dir + 'Modeling/cleaned_template_fastas/'):
        os.makedirs(parent_dir + 'Modeling/cleaned_template_fastas/')
    if not os.path.exists(parent_dir + 'Modeling/cleaned_template_pdbs/'):
        os.makedirs(parent_dir + 'Modeling/cleaned_template_pdbs/')
    if not os.path.exists(parent_dir + 'Modeling/fasta_alns_and_identities/'):
        os.makedirs(parent_dir + 'Modeling/fasta_alns_and_identities/')
    if not os.path.exists(parent_dir + 'Modeling/grishin_alns/'):
        os.makedirs(parent_dir + 'Modeling/grishin_alns/')
    if not os.path.exists(parent_dir + 'Modeling/threaded_pdbs/'):
        os.makedirs(parent_dir + 'Modeling/threaded_pdbs/')
    if not os.path.exists(parent_dir + 'Modeling/final_models/'):
        os.makedirs(parent_dir + 'Modeling/final_models/')


def get_chains_from_swiss_model(path_to_swiss_model):
    """
    Parses chains in pdb file with swiss model

    :param path_to_swiss_model: path to pdb file (ex: path/to/file.pdb)
    :type path_to_swiss_model: str
    :return: letters corresponding to all chains found in pdb file
    :rtype: list
    """

    with open(path_to_swiss_model) as pdb:
        line = 'empty'
        uniq_chains = []
        while line:
            line = pdb.readline()
            if line[0:4] == 'ATOM':
                amino_acids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN',
                               'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
                line_contents = line.split(' ')
                for i in range(len(line_contents)):
                    if line_contents[i] in amino_acids:
                        chain = line_contents[i + 1][0]
                        if chain not in uniq_chains:
                            uniq_chains.append(chain)
    return uniq_chains


def clean_pdb_output_of_clean_pdb(path, template, chain):
    """
    As part of the run_clean_pdb, this function takes its output pdb file as input (located in provided directory,
    contains one chain template), corrects its name if necessary and moves it to cleaned_template_pdbs directory \n
    Returns nothing

    :param path: path to pdb file (ex: path/to/file.pdb)
    :type path: str
    :param template: pdb file current name
    :type template: str
    :param chain: the selected chain of a template
    :type chain: str
    """

    if template[0].isdigit():  # crystal
        new_pdb = template[:-4] + '.pdb'
        os.rename(template[:-4] + '_' + chain + '.pdb', new_pdb)  # one chain in name is enough
        shutil.move(new_pdb, path + 'Modeling/cleaned_template_pdbs')
    else:  # swiss_model
        shutil.move(template[:-4] + '_' + chain + '.pdb', path + 'Modeling/cleaned_template_pdbs')


def clean_fasta_output_of_clean_pdb(path, template, chain):
    """
    As part of the run_clean_pdb, this function takes its output fasta file as input (located in provided directory,
    contains one chain template), corrects comment line in file, file's name if necessary and moves it to
    cleaned_template_fasts directory \n
    Returns nothing

    :param path: path to fasta file (ex: path/to/file.fasta)
    :type path: str
    :param template: fasta file current name
    :type template: str
    :param chain: the selected chain of a template
    :type chain: str
    """

    if template[0].isdigit():                             # template is from wwPDB db
        new_fasta = template[:-4] + '.fasta'
        os.rename(template[:-4] + '_' + chain + '.fasta', new_fasta)
        # correct name inside file
        with open(new_fasta, 'r') as file:
            content = file.readlines()
            content[0] = content[0][:-3]
        os.remove(new_fasta)
        with open(path + 'Modeling/cleaned_template_fastas/' + new_fasta, 'w') as file:
            for line in content:
                file.write(line)
                if '>' in line:
                    file.write('\n')
    else:                                                 # template is from Swiss-model db
        shutil.move(template[:-4] + '_' + chain + '.fasta', path + 'Modeling/cleaned_template_fastas/')


def run_clean_pdb(working_dir, path_to_rosetta_script):
    """
    Executes Rosetta' clean_pdb.py script on all selected templates. Differentiates between structures from Swiss-Model
    and wwPDB databases (unlike Swiss-model structures wwPDB structure files in this script contain the selected chain
    letter in its name as last symbol before .pdb format) \n
    Returns nothing

    :param working_dir: path to fasta file (ex: path/to/directory/)
    :type working_dir: str
    :param path_to_rosetta_script: path to directory with compiled Rosetta (ex: path/to/Rosetta/)
    :type path_to_rosetta_script: str
    """

    path_to_templates = working_dir + 'orig_templates/'
    templates = next(os.walk(path_to_templates))[2]
    print(f'\nFound templates: {templates}\n')
    for i in range(len(templates)):
        if templates[i][0].isdigit():
            print(f'\nTemplate {templates[i]} is from wwPDB database')
            chain = templates[i][-5]  # chain is the last symbol of crystal file name ([-5] to ignore .pdb)
            print(f'The chosen chain is {chain}')
            subprocess.run(
                ['python3', path_to_rosetta_script + 'clean_pdb.py', path_to_templates + templates[i], chain])
            print('Cleaning pdb file')
            clean_pdb_output_of_clean_pdb(working_dir, templates[i], chain)
            print('Cleaning fasta file')
            clean_fasta_output_of_clean_pdb(working_dir, templates[i], chain)
        else:
            print(f'\nTemplate {templates[i]} is from Swiss-model database\n')
            chains = get_chains_from_swiss_model(path_to_templates + templates[i])
            print(f'\nThe template contains chains {chains}\n')
            for chain in chains:
                print(f'\nRun clean_pdb.py on {templates[i]} with chain {chain}\n')
                subprocess.run(
                    ['python3', path_to_rosetta_script + 'clean_pdb.py', path_to_templates + templates[i], chain])
                print('Cleaning pdb file\n')
                clean_pdb_output_of_clean_pdb(working_dir, templates[i], chain)
                print('\nCleaning fasta file\n')
                clean_fasta_output_of_clean_pdb(working_dir, templates[i], chain)


def mafft_multiple_alignment(path, id_protein, output_name):
    """
    Generates fasta alignment by running mafft L-INS-i alignment on target protein (id_protein) with fasta files
    gathered from the provided directory \n
    Returns number of aligned sequences

    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param id_protein: Uniprot ID of target protein (ex: Q86WT6)
    :type id_protein: str
    :param output_name: output fasta alignment file name
    :type output_name: str
    :return: number of aligned fasta sequences
    :rtype: int
    """

    path_to_templates = path + 'Modeling/cleaned_template_fastas/'
    path_to_target = path + id_protein + '.fasta'
    with open('fastas_for_mafft', 'w') as fastas:

        # write target fasta in joint file

        target = open(path_to_target)
        for line in target:
            fastas.write(line)
        fastas.write(line)
        target.close()

        # write templates fastas in joint file

        number_of_fastas = 1  # 1 is for target
        templates = next(os.walk(path_to_templates))[2]
        print(templates)
        for i in templates:
            number_of_fastas += 1
            with open(path_to_templates + i) as template:
                for line in template:
                    fastas.write(line)
    path_to_alignment = path + 'Modeling/fasta_alns_and_identities/'
    os.system('mafft --localpair --maxiterate 1000 fastas_for_mafft > ' + path_to_alignment + output_name)
    # os.remove('fastas_for_mafft')
    return number_of_fastas


def remove_duplicates(file, number_of_fastas, path, output_name):
    """
    Takes fasta alignment file and generates its copy free of identical fasta sequences from original alignment and
    deletes corresponding pdb and fasta files from cleaned_template_pdbs and cleaned_template_fastas directories,
    respectively. Moves the resulted alignment to fasta_alns_and_identities directory \n
    Returns number of unique fasta sequences

    :param file: file with multiple fasta sequences
    :type file: str
    :param number_of_fastas: number of fasta sequences in input file
    :type number_of_fastas: int
    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param output_name: output fasta alignment file name
    :type output_name: str
    :return: number of unique fasta sequences
    :rtype: int
    """

    path_to_pbds = path + 'Modeling/cleaned_template_pdbs/'
    path_to_fastas = path + 'Modeling/cleaned_template_fastas/'
    path_to_alignnment = path + 'Modeling/fasta_alns_and_identities/' + file
    fastas = parse_multifasta_file(path_to_alignnment, number_of_fastas)
    uniq_fastas = []
    with open(output_name, "w") as f:
        for i in range(number_of_fastas):
            name, seq = next(fastas)
            if seq not in uniq_fastas:
                uniq_fastas.append(seq)
                f.write('>' + name + '\n')
                f.write(seq + '\n')
            else:
                os.remove(path_to_pbds + name + '.pdb')
                os.remove(path_to_fastas + name + '.fasta')
    shutil.move(output_name, path + 'Modeling/fasta_alns_and_identities/')
    return len(uniq_fastas)


def calculate_identity(target, template):
    """
    Calculates identity percent between aligned target and template sequences (ignores cases when both sequences
    contain '-' symbol)

    :param target: target sequence
    :type target: str
    :param template: template sequence
    :type template: str
    :return: identity percent
    :rtype: float
    """

    identical_residues = 0
    length_of_target = 0
    for i in range(len(target)):
        if target[i] != '-':
            length_of_target += 1
            if target[i] == template[i]:
                identical_residues += 1
    identity_percent = round(identical_residues / length_of_target * 100, 2)
    return identity_percent


def calculate_identity_for_mult_templates(path, alignment, number_of_fastas):
    """
    Calculates identity percent between target sequence (first fasta in file) and every aligned template sequence
    from alignment file. Generates file with calculated values in fasta_alns_and_identities directory \n
    Returns nothing

    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param alignment: alignment file
    :type alignment: str
    :param number_of_fastas: number of fasta sequences in input file
    :type number_of_fastas: int
    """

    path_to_alignment = path + 'Modeling/fasta_alns_and_identities/' + alignment
    fastas = parse_multifasta_file(path_to_alignment, number_of_fastas)
    target_name, target_seq = next(fastas)
    with open('identity', 'w') as iden_file:
        for i in range(number_of_fastas - 1):
            homolog_name, homolog_seq = next(fastas)
            identity = calculate_identity(target_seq, homolog_seq)
            iden_file.write(homolog_name + ' ' + str(identity) + '\n')
    shutil.move('identity', path + 'Modeling/fasta_alns_and_identities/')


def remove_templates_under_threshold(path, identity_file, threshold):
    """
    Takes identity percents generated by calculate_identity_for_mult_templates function and removes all templates in
    cleaned_template_pdbs and cleaned_template_fastas directories with identity percent below threshold \n
    Returns nothing

    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param identity_file: file with identity percents
    :type identity_file: str
    :param threshold: minimum identity percent
    :type threshold: float
    """

    path_to_pbds = path + 'Modeling/cleaned_template_pdbs/'
    path_to_fastas = path + 'Modeling/cleaned_template_fastas/'
    path_to_alignment = path + 'Modeling/fasta_alns_and_identities/' + identity_file
    to_remove = []
    with open(path_to_alignment) as identities:
        identity = 'empty'
        while identity[0]:
            identity = identities.readline().strip().split(' ')
            if identity[0]:
                if float(identity[1]) < threshold:
                    to_remove.append(identity[0])
    for i in to_remove:
        os.remove(path_to_pbds + i + '.pdb')
        os.remove(path_to_fastas + i + '.fasta')


def calculate_coverage(path, alignment, number_of_fastas):
    """
    Calculates coverage percent - percent of target sequence covered by all selected aligned templates together \n
    Returns coverage percent

    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param alignment: alignment file
    :type alignment: str
    :param number_of_fastas: number of fasta sequences in input file
    :type number_of_fastas: int
    :return: coverage percent
    :rtype: float
    """

    path_to_alignment = path + 'Modeling/fasta_alns_and_identities/' + alignment
    fastas_iterator = parse_multifasta_file(path_to_alignment, number_of_fastas)
    fastas = []
    targer_name, target_seq = next(fastas_iterator)
    fastas.append(target_seq)
    length_of_target = 0
    for i in target_seq:
        if i != '-':
            length_of_target += 1
    for i in range(1, number_of_fastas):
        name, seq = next(fastas_iterator)
        fastas.append(seq)
    coverage = 0
    for i in range(len(fastas[0])):
        for j in range(1, len(fastas)):
            if fastas[0][i] != '-' and fastas[j][i] != '-':
                coverage += 1
                break
    coverage_percent = round(coverage / length_of_target * 100, 2)
    return coverage_percent


def convert_fasta_alingment_to_grishin(path, alignment, number_of_fastas):
    """
    Converts multiple fasta alignment into paired Grishin alignments required for threading and moves generated files
    to grishin_alns directory \n
    Returns nothing

    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param alignment: alignment file
    :type alignment: str
    :param number_of_fastas: number of fasta sequences in input file
    :type number_of_fastas: int
    """

    path_to_alignment = path + 'Modeling/fasta_alns_and_identities/' + alignment
    fastas = parse_multifasta_file(path_to_alignment, number_of_fastas)
    target_name, target_seq = next(fastas)
    for i in range(1, number_of_fastas):
        homolog_name, homolog_seq = next(fastas)
        with open(target_name + '_' + homolog_name + '.grishin', 'w') as grishin_aln:
            (grishin_aln.write('## ' + target_name + ' ' + homolog_name + '.pdb\n#\nscores from program: 0\n0 ' +
                               target_seq + '\n0 ' + homolog_seq))
            grishin_aln.close()
            shutil.move(target_name + '_' + homolog_name + '.grishin', path + 'Modeling/grishin_alns')


def thread_selected_pdbs(path, target, path_to_rosetta):
    """
    Runs Rosetta partial thread application on target sequence over all selected templates individually. Moves generated
    pdb file to threaded_pdbs directory \n
    Returns nothing

    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param target: Uniprot ID of target protein (ex: Q86WT6)
    :type target: str
    :param path_to_rosetta: path to directory with compiled Rosetta (ex: path/to/Rosetta/)
    :type path_to_rosetta: str
    """

    path_to_target = path + target + '.fasta'
    path_to_grishin_alins = path + 'Modeling/grishin_alns/'
    path_to_template_pdbs = path + 'Modeling/cleaned_template_pdbs/'
    path_to_rosetta_script = path_to_rosetta + 'main/source/bin/partial_thread.default.linuxgccrelease'
    templates = next(os.walk(path_to_template_pdbs))[2]
    grishins = next(os.walk(path_to_grishin_alins))[2]
    for template in templates:
        template_name = template[:-4]
        template_r = re.compile('.*({}).*'.format(template_name))
        grishin = list(filter(template_r.match, grishins))[0]
        print(grishin)
        subprocess.run([path_to_rosetta_script,
                        '-in:file:fasta ' + path_to_target,
                        '-in:file:alignment ' + path_to_grishin_alins + grishin,
                        '-in:file:template_pdb ' + path_to_template_pdbs + template])
        new_threaded_file = (template[:-4] + '_on_' + target + '.pdb')
        os.rename(template[:-4] + '.pdb.pdb', new_threaded_file)
        shutil.move(new_threaded_file, path + 'Modeling/threaded_pdbs/')


def make_fragments(path, target, path_to_rosetta):
    """
    Runs Rosetta make_fragments.pl script on target sequence. Moves output fragment file to Modeling directory \n
    Returns nothing

    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param target: Uniprot ID of target protein (ex: Q86WT6)
    :type target: str
    :param path_to_rosetta: path to directory with compiled Rosetta (ex: path/to/Rosetta/)
    :type path_to_rosetta: str
    """

    path_to_target = path + target + '.fasta'
    path_to_rosetta_script = path_to_rosetta + 'main/tools/fragment_tools/make_fragments.pl'
    target = parse_multifasta_file(path_to_target, 1)
    target_name, target_seq = next(target)
    subprocess.run([path_to_rosetta_script, '-verbose', path_to_target])
    os.rename('aat000_03_05.200_v1_3', target_name + '_3.frags')
    os.rename('aat000_09_05.200_v1_3', target_name + '_9.frags')
    shutil.move(target_name + '_3.frags', path + 'Modeling/')
    shutil.move(target_name + '_9.frags', path + 'Modeling/')


def modify_xml_file(path, target):
    """
    Modifies hybridize_empty.xml file (protocol for Rosetta cm) by adding corresponding fragment and threaded pdb files,
    moves it to Modeling directory

    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param target: Uniprot ID of target protein (ex: Q86WT6)
    :type target: str
    """

    path_to_template_pdbs = path + 'Modeling/threaded_pdbs/'
    templates = next(os.walk(path_to_template_pdbs))[2]
    parser = etree.XMLParser(remove_blank_text=True)
    tree = etree.parse('hybridize_empty.xml', parser)
    root = tree.getroot()
    hybridize = root[3][0]
    fragments = root[3][0][0]

    three_mers_name = path + 'Modeling/' + target + '_3.frags'
    nine_mers_name = path + 'Modeling/' + target + '_9.frags'
    fragments.set("three_mers", three_mers_name)
    fragments.set("nine_mers", nine_mers_name)

    for i in range(len(templates)):
        template = etree.Element("Template")
        template.set("pdb", path_to_template_pdbs + templates[i])
        template.set('cst_file', "AUTO")
        template.set('weight', '1.000')
        hybridize.append(template)

    tree.write("hybridize.xml", pretty_print=True)
    shutil.move("hybridize.xml", path + 'Modeling/')


def run_modeling(path, path_to_rosetta, target, nsrtuct):
    """
    Runs The Multi-template Comparative Rosetta Modeling. The resulted models are renamed to contain target ID and
    moved to final_models directory \n
    Returns nothing

    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param path_to_rosetta: path to directory with compiled Rosetta (ex: path/to/Rosetta/)
    :type path_to_rosetta: str
    :param target: Uniprot ID of target protein (ex: Q86WT6)
    :type target: str
    :param nsrtuct: the desired number of final models
    :type nsrtuct: str
    """

    path_to_rosetta_script = path_to_rosetta + 'main/source/bin/rosetta_scripts.default.linuxgccrelease'
    path_to_rosetta_database = path_to_rosetta + 'main/database'
    path_to_target = path + target + '.fasta'
    path_to_protocol = path + 'Modeling/hybridize.xml'
    subprocess.run([path_to_rosetta_script,
                    '-database', path_to_rosetta_database,
                    '-in:file:fasta', path_to_target,
                    '-parser:protocol', path_to_protocol,
                    '-default_max_cycles', '200',
                    '-dualspace',
                    '-nstruct', nsrtuct,
                    '-restore_talaris_behavior',
                    '-score:set_weights', 'pro_close', '0', 'cart_bonded', '0.5'])
    for i in range(1, int(nsrtuct) + 1):
        name = 'modeled_' + target + '_' + str(i) + '.pdb'
        os.rename('S_000' + str(i) + '.pdb', name)
        shutil.move(name, path + 'Modeling/final_models/')


