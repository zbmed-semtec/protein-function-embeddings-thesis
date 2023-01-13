import time
import pandas as pd
import numpy as np
import multiprocessing


def read_embeddings(input_embeds_file: str) -> np.array:
    """
    Reads the embeddings pickle file and returns a numpy array of all accession ids.
    Parameters
    ----------
    input_embeds_file: str
        Filepath to the embeddings pickle file (accepts pickle file for any hyperparameter combination)
    Returns
    -------
    all_accessions : np.array
        Numpy array consisting of all accessions.
    """
    data = pd.read_pickle(input_embeds_file)
    all_accessions = np.array(data['accessions'].tolist())
    return all_accessions


def read_uniref_clusters(input_cluster_file: str) -> np.array:
    """
    Reads the Uniref cluster TSV file and returns a numpy array of cluster members.
    Parameters
    ----------
    input_cluster_file : str
        Filepath to Uniref clusters TSV file.
    Returns
    -------
    uniref_cluster_members : np.array
        Nested numpy array of cluster members.
    """
    data = pd.read_csv(input_cluster_file, sep='\t', names=['cluster_id', 'cluster_name', 'members', 'member_count'])
    to_array = lambda x: np.array(eval(x))
    data['members'] = data['members'].apply(to_array)
    uniref_cluster_members = np.array(data['members'])
    return uniref_cluster_members


def extract_accession(x: int, all_accessions: np.array) -> str:
    """
    Returns the accession number at index 'x' of the input all_accessions array.
    Parameters
    ----------
    x : int
        Index at which the accession is to be retrieved.
    all_accessions : np.array
        Numpy array consisting of all accessions.
    Returns
    -------
    all_accessions[x] : str
        Accession at index 'x'.
    """
    return all_accessions[x]


def calculate_true_positive(pair: list) -> list[int, str, int, list]:
    """
    Check if the accession belongs to any Uniref clusters. Returns the index of accession, accession, whether it is
    present or not and index of cluster where/if it is present.
    Parameters
    ----------
    pair : list
        List of accession index and accession.
    Returns
    -------
    accession_index : int
        Index of accession.
    accession : str
        Uniprot accession id.
    presence : int
        If accession belong to any cluster or not (1 = present, 0 = not present)
    index : List[int]
        List of cluster indices where the accession is present.
    """
    accession_index = pair[0]
    accession = pair[1]
    values = [accession]
    result = [len(set(values).intersection(uniref_clusters[i])) for i in range(uniref_clusters.shape[0])]
    if result.count(1) > 0:
        presence = 1
        index = [i for i, x in enumerate(result) if x == 1]
    else:
        presence = 0
        index = [0]
    return [accession_index, accession, presence, index]


if __name__ == "__main__":

    all_accessions = read_embeddings("./data/word2doc2vec_embs_euk_1.pkl")
    print("Read accessions of embeddings")

    uniref_clusters = read_uniref_clusters("./data/uniref_clusters.tsv")
    print(f"Read clusters, total clusters: {len(uniref_clusters)}")

    accessions = []
    for idx, accession_index in enumerate(all_accessions):
        accession = extract_accession(idx, all_accessions)
        accessions.append([idx, accession])
    print(f"Total pairs {len(accessions)}")

    print("Starting index fetch")
    start_time = time.time()
    with multiprocessing.Pool(10) as pool:
        results = pool.map(calculate_true_positive, accessions)
    df = pd.DataFrame(results)
    df.to_csv("./data/cluster_index.tsv", sep='\t', index=False, header=False)
    print("--- %s seconds mp---" % (time.time() - start_time))


