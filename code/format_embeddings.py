import os
import numpy as np
import pandas as pd
from typing import Tuple


def load_embeddings(embeddings_dir: str) -> Tuple[list, list]:
    """
    Loads each NPY file in the given embeddings directory and returns a list of all accession ids and embeddings.
    Parameters
    ----------
    embeddings_dir : str
        Directory path to the embeddings.
    Returns
    -------
    accessions : list
        List of all accession ids.
    embeddings : list
        List of all embeddings.
    """
    accessions = []
    embeddings = []
    for embedding_file in os.listdir(embeddings_dir):
        accession = embedding_file.split(".npy")[0]
        embedding = np.load(f'{embeddings_dir}/{embedding_file}', allow_pickle=True)
        accessions.append(accession)
        embeddings.append(np.array(embedding, dtype=np.float16))
    return accessions, embeddings


def save_embeddings(accessions: list, embeddings: list, output_pickle_path: str) -> None:
    """
    Saves the embeddings in the pickle format with sorted accession ids as the first column and embeddings as the second
    column.
    Parameters
    ----------
    accessions : list
        List of all accession ids.
    embeddings : list
        List of all embeddings.
    output_pickle_path : str
        File path for the output pickle file.
    """
    embeddings_df = pd.DataFrame({"accessions": accessions, "embeddings": embeddings})
    embeddings_df.sort_values(by="accessions", ignore_index=True, inplace=True)
    embeddings_df.to_pickle(output_pickle_path)


def load_embeddings_pickle(input_pickle_path: str) -> pd.DataFrame:
    """
    Loads and returns a dataframe from the embeddings pickle file.
    Parameters
    ----------
    input_pickle_path

    Returns
    -------
    unpickled_embeddings_df : pd.DataFrame
        Pandas Dataframe of accession ids and embeddings.
    """
    unpickled_embeddings_df = pd.read_pickle(input_pickle_path)
    return unpickled_embeddings_df


if __name__ == "__main__":
    all_accessions, all_embeddings = load_embeddings("./data/output/embeddings/word2doc2vec/cbow/min_count_2")
    save_embeddings(all_accessions, all_embeddings, "./data/output/embeddings_pickle/word2doc2vec_embs_1.pkl")
    load_embeddings_pickle("./data/output/embeddings_pickle/word2doc2vec_embs_1.pkl")