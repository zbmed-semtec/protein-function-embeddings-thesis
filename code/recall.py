import warnings
import csv
import pandas as pd
import numpy as np
from tqdm import tqdm
from numba_progress import ProgressBar


def read_accessions(embeds_path: str) -> np.array:
    """
    Reads the embeddings pickle file as a pandas dataframe and returns the accessions column as a numpy array.
    Parameters
    ----------
    embeds_path : str
        Filepath to the embeddings pickle file.
    Returns
    -------
    all_accessions : np.array
        Numpy array of all Eukaryota accessions.
    """
    data = pd.read_pickle(embeds_path)
    all_accessions = np.array(data['accessions'].tolist())
    return all_accessions


def read_test_accessions(test_accession_path: str) -> np.array:
    """
    Reads the test accessions (20% of Eukaryota entries) as a pandas dataframe and returns the index column as a numpy
     array. The index column is the corresponding index value based on all Eukaryota entries.
    Parameters
    ----------
    test_accession_path : str
        Filepath for 20% of Eukaryota entries.
    Returns
    -------
    test_accessions_reference : np.array
        Numpy array of index column of test accessions (20% Eukaryota entries).
    """
    test_data = pd.read_csv(test_accession_path, sep='\t')
    test_accessions_reference = np.array(test_data["index"].tolist())
    return test_accessions_reference


def read_cluster_index(cluster_ids: str) -> set:
    """
    Reads the file containing clustered pairs along with the index values as a pandas dataframe and returns a
    set of tuple containing the index values of both the accessions of a pair.
    Parameters
    ----------
    cluster_ids : str
        Filepath to the clustered pairs along with the indices.
    Returns
    -------
    ids: set
        Set of tuple of index value of accession1 and accession2 of every pair.
    """
    all_ids = pd.read_csv(cluster_ids, sep='\t')
    all_ids = all_ids.apply(lambda row: [row['accession1_index'], row['accession2_index']], axis=1).to_numpy()
    ids = set(map(tuple, all_ids))
    return ids


def read_cosine_matrix(matrix_file: str) -> np.matrix:
    """
    Reads, loads and returns the cosine matrix.
    Parameters
    ----------
    matrix_file : str
        Filepath to the cosine matrix.
    Returns
    -------
    cosine_matrix : np.matrix
        Numpy matrix of the cosine matrix.
    """
    matrix = np.load(matrix_file)
    cosine_matrix = matrix['arr_0']
    return cosine_matrix


def extract_test_accession_index(test_accession_index: int, test_accessions_reference: np.array) -> int:
    """
    Returns the index value of the test accession based on all Eukaryota entries.
    Parameters
    ----------
    test_accession_index : int
        Index value of test accession based on 20% of the entries.
    test_accessions_reference : np.array
        Numpy array of index column of test accessions (20% Eukaryota entries).
    Returns
    -------
    test_accessions_reference[test_accession_index] : int
        Index value of test accession based on its position in all Eukaryota entries.
    """
    return test_accessions_reference[test_accession_index]


def calculate_recall(cosine_matrix: np.matrix, progress_proxy: tqdm, subset: bool = True) -> tuple[int, int]:
    """
    Parses the cosine matrix, extracts all pairs having a cosine score of >= 0.90 and checks if the pair belongs to a
    clustered pair. If the pair belongs to a clustered pair, it is counted as a True positive.
    Returns the count of all true positives.
    Parameters
    ----------
    cosine_matrix : np.matrix
        Numpy cosine matrix.
    progress_proxy : tqdm
        tdqm progress bar.
    subset : bool
        If true, calculates recall for subset (20% of Eukaryota entries). Otherwise, for all Eukaryota entries.
    Returns
    -------
    true_positives : int
        Count of true positives.
    """
    rows = (cosine_matrix != 0).sum()
    print(f"Total 0.90 pairs {rows}")
    true_positives = 0
    for test_accession_index, row in enumerate(cosine_matrix):
        cosine_90_indices = np.nonzero(row)[0]
        for idx, accession_index in enumerate(cosine_90_indices):
            if subset:
                test_ref_index = extract_test_accession_index(test_accession_index, test_accessions_reference)
                pair1 = tuple(np.array([accession_index, test_ref_index]))
                pair2 = tuple(np.array([test_ref_index, accession_index]))
            else:
                pair1 = tuple(np.array([accession_index, test_accession_index]))
                pair2 = tuple(np.array([test_accession_index, accession_index]))
            if pair1 in ids or pair2 in ids:
                    true_positives += 1
        progress_proxy.update(1)
    false_positives = len(ids) - true_positives
    return true_positives, false_positives


def write_recall_scores(hyperparameter: str, true_positives: int, false_positives: int, recall_file: str) -> None:
    """

    Parameters
    ----------
    hyperparameter : str
        Name of hyperparameter combination.
    true_positives : int
        Count of true positives.
    false_positives : int
        Count of false negatives.
    recall_file : str
        Filepath to store recall scores.
    """
    headers = ['Parameter', "True positives", "False positves", "Recall"]
    data = [hyperparameter, true_positives, false_positives, true_positives/(true_positives+false_positives)]
    with open(recall_file, 'w', newline='') as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerow(headers)
        writer.writerow(data)


if __name__ == "__main__":
    warnings.filterwarnings('ignore')

    embeds_path = "./data/emb_pickle/word2doc2vec_embs_euk_1.pkl"
    all_accessions = read_accessions(embeds_path)
    print("Read accessions")

    test_accession_path = "./data/files/rev-20220525-UniProtKB-eukaryota-20.tsv"
    test_accessions_reference = read_test_accessions(test_accession_path)
    print("Read reference accessions")

    cluster_ids = "./data/clustered_pairs_index.tsv"
    ids = read_cluster_index(cluster_ids)
    print("Read ids")

    matrix_file = "./data/cosine/cosine_word2doc2vec_1.npz"
    cosine_matrix = read_cosine_matrix(matrix_file)
    print("Loaded cosine matrix")

    with ProgressBar(total=cosine_matrix.shape[0]) as progress:
        true_positives, false_positives = calculate_recall(cosine_matrix, progress, subset=False)
        write_recall_scores("word2doc2vec_1", true_positives, false_positives, "./data/recall_scores.tsv")
    print("Calculated recall")

