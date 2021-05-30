import shutil
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os


def make_output_dirs_for_part3(parent_dir):
    """
    Functions creates score directory necessary for the third part 'Scoring' in provided working directory
    (takes into account the case if the directory exists already) \n
    Returns nothing

    :param parent_dir: path to working directory (ex: path/to/directory/)
    :type parent_dir: str
    """

    if not os.path.exists(parent_dir + 'score/'):
        os.makedirs(parent_dir + 'score/')


def run_ornate_score(path, target, path_to_ornate):
    """
    Runs Ornate scoring for Rosetta models. Generated score files are moved to score directory \n
    Returns nothing

    :param path: path to working directory (ex: path/to/directory/)
    :type path: str
    :param target: Uniprot ID of target protein (ex: Q86WT6)
    :type target: str
    :param path_to_ornate: path to directory with Ornate tool (ex: path/to/Ornate/)
    :type path_to_ornate: str
    """

    path_to_models = path + 'Modeling/final_models/'
    models = next(os.walk(path_to_models))[2]
    for i in range(1, len(models) + 1):
        score_file_name = 'score_modeled_' + target + '_' + str(i) + '.txt'
        os.chdir(path_to_ornate)
        os.system(f'python3.7 score.py -s {path_to_models + models[i - 1]} > {score_file_name}')
        os.chdir(path)
        shutil.move(path_to_ornate + score_file_name, path + 'score/')


def easy_parser(name):
    '''
    Parse ".txt" files with score \n
    Returns pd.DataFrame

    :param name: ".txt" file
    :type name: str
    :return:pd.DataFrame
    '''
    colnames=['.', 'number', 'res', 'score']
    df = pd.read_csv(f"{name}", sep = "\s+", names=colnames, header=None)
    df = df.drop([0, 1], axis='index')
    df = df.drop(".", axis = 1)
    df["score"] = pd.to_numeric(df["score"])
    return df


def score_counter(directory, mean_and_sd=True, graph=True):
    '''
    Counter mean and sd score value for each model \n
    Returns score files ".csv" and/or ".jpg"".

    :param directory: dir where the db file will be saved
    :type directory: str
    :param mean_and_sd: if choose this flag function will be count mean and sd, and return ".csv"
    :type direcmean_and_sd: bool
    :param graph:if choose this flag function will be to create plot with the a quality score for the residuals".jpg"
    :return: files ".csv" and/or ".jpg"
    :rtype:file
    '''
    path_to_score = os.path.join(directory, "score")
    os.chdir(path_to_score)
    text_files = [f for f in os.listdir(path_to_score) if f.endswith('.txt')]
    name = []
    mean = []
    sd = []
    fig = plt.figure(num=None, figsize=(10, 4), dpi=100, facecolor='w', edgecolor='k')
    ax = plt.axes()
    ax.grid()
    ax.set_xlabel('Number_of_residual')
    ax.set_ylabel('Score')
    count = 0
    for i in text_files:
        name.append(i[:-4])
        df_model = easy_parser(i)
        mean.append(round(df_model["score"].mean(), 3))
        sd.append(round(df_model["score"].std(), 3))
        plt.plot(df_model["number"], df_model["score"], label = name[count])
        count +=1
        length_sequence = len(i)
    score_result = os.path.join(directory, "summary_result")
    score_result_to_csv = os.path.join(score_result, "score_mean_sd.csv")
    score_result_graph = os.path.join(score_result, "graph_with_res.jpg")
    if not os.path.exists(score_result):
        os.makedirs(score_result)
    if mean_and_sd:
        df_metrics = pd.DataFrame({'name': name,
                   'mean': mean,
                   'sd': sd})
        df_metrics.to_csv(score_result_to_csv, index=False)
    step = 25
    plt.xticks(np.arange(0, length_sequence, step))
    plt.legend()
    if graph:
        plt.savefig(score_result_graph)
    print("Success !")