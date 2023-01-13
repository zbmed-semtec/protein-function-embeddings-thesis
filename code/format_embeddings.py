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
        embeddings.append(np.array(embedding, dtype=np.float32))
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
    embeddings_df.drop(embeddings_df[embeddings_df['accessions'] == 'P19664'].index, inplace=True)
    embeddings_df.drop(embeddings_df[embeddings_df['accessions'] == 'P22972'].index, inplace=True)
    embeddings_df.to_pickle(output_pickle_path, protocol=4)
    print("Saved embeddings")


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

    # WORD2DOC2VEC 200
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/200/cbow/min_count_2/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_1.pkl")
    print("Embeddings formatted 1")
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/200/cbow/min_count_3/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_2.pkl")
    print("Embeddings formatted 2")
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/200/cbow/min_count_4/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_3.pkl")
    print("Embeddings formatted 3")
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/200/sg/min_count_2/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_4.pkl")
    print("Embeddings formatted 4")
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/200/sg/min_count_3/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_5.pkl")
    print("Embeddings formatted 5")
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/200/sg/min_count_4/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_6.pkl")
    print("Embeddings formatted 6")

    # WORD2DOC2VEC 400
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/400/cbow/min_count_2/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_7.pkl")
    print("Embeddings formatted 7")
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/400/cbow/min_count_3/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_8.pkl")
    print("Embeddings formatted 8")
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/400/cbow/min_count_4/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_9.pkl")
    print("Embeddings formatted 9")
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/400/sg/min_count_2/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_10.pkl")
    print("Embeddings formatted 10")
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/400/sg/min_count_3/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_11.pkl")
    print("Embeddings formatted 11")
    all_accessions, all_embeddings = load_embeddings("./data/embeddings/word2doc2vec/400/sg/min_count_4/")
    save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/word2doc2vec_embs_euk_12.pkl")
    print("Embeddings formatted 12")

    # HYBRID 200
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/200/cbow/min_count_2/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_1.pkl")
    # print("Embeddings formatted 1")
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/200/cbow/min_count_3/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_2.pkl")
    # print("Embeddings formatted 2")
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/200/cbow/min_count_4/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_3.pkl")
    # print("Embeddings formatted 3")
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/200/sg/min_count_2/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_4.pkl")
    # print("Embeddings formatted 4")
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/200/sg/min_count_3/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_5.pkl")
    # print("Embeddings formatted 5")
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/200/sg/min_count_4/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_6.pkl")
    # print("Embeddings formatted 6")

    # HYBRID 400
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/400/cbow/min_count_2/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_7.pkl")
    # print("Embeddings formatted 7")
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/400/cbow/min_count_3/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_8.pkl")
    # print("Embeddings formatted 8")
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/400/cbow/min_count_4/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_9.pkl")
    # print("Embeddings formatted 9")
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/400/sg/min_count_2/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_10.pkl")
    # print("Embeddings formatted 10")
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/400/sg/min_count_3/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_11.pkl")
    # print("Embeddings formatted 11")
    # all_accessions, all_embeddings = load_embeddings("./data/embeddings/hybrid/400/sg/min_count_4/")
    # save_embeddings(all_accessions, all_embeddings, "./data/emb_pickle/hybrid_embs_euk_12.pkl")
    # print("Embeddings formatted 12")