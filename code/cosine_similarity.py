import numpy as np
import pandas as pd
from numba import njit
from tqdm import tqdm
from numba_progress import ProgressBar


@njit(fastmath=True)
def cosine_similarity_numba(u: np.ndarray, v: np.ndarray) -> np.float64:
    """
    Computes the cosine similarity between two vectors using numba
    Parameters
    ----------
    u : np.ndarray
        First/Reference embedding.
    v : np.ndarray
        Second/Assesed embedding.
    Returns
    -------
    cos_theta : np.float64
        Cosine similarity score between two embeddings.
    """
    assert(u.shape[0] == v.shape[0])

    uv = 0
    uu = 0
    vv = 0
    for i in range(u.shape[0]):
        uv += u[i]*v[i]
        uu += u[i]*u[i]
        vv += v[i]*v[i]
    cos_theta = 1
    if uu != 0 and vv != 0:
        cos_theta = uv/np.sqrt(uu*vv)
    return cos_theta


@njit(nogil=True)
def compute_matrix(embeds: np.array, test_data_indices: np.array, progress_proxy: tqdm) -> np.array:
    """
    Creates an empty numpy matrix with rows based on the dimensions of input test_data_indices and columns based on the
    dimensions of embeddings. Iterates through each cell of the matrix to calculate the cosine similarity between the
    pair of embeddings. Stores the cosine similarity score if it is greater or equal to 0.90 and returns the filled
    matrix.
    Parameters
    ----------
    embeds : np.array
        All SwissProt Eukaryota embeddings in a numpy array.
    test_data_indices : np.array
        All indices of the 20 percent of SwissProt Eukaryota embeddings in a numpy array.
    progress_proxy : tqdm
        tqdm progress bar proxy object to show progress.
    Returns
    -------
    matrix : Cosine similarity matrix.
    """
    matrix = np.zeros((test_data_indices.shape[0], embeds.shape[0]), dtype=np.float32)
    for i, index in enumerate(test_data_indices):
        embed1 = embeds[index]
        for j in range(embeds.shape[0]):
            embed2 = embeds[j]
            cs = round(cosine_similarity_numba(embed1, embed2), 4)
            if cs >= 0.9:
                matrix[i, j] = cs
        progress_proxy.update(1)
    return matrix


def create_matrix(embeds_path: str, test_accession_path: str, matrix_filepath: str) -> None:
    """
    Wrapper function to read the SwissProt Eukaryota embeddings pickle file and 20 percent of SwissProt Eukarota entries file,
    compute cosine similarity and save the cosine similarity matrix as a .NPZ file.
    Parameters
    ----------
    embeds_path : str
        Filepath to SwissProt Eukaryota embeddings pickle file.
    test_accession_path : str
        Filepath to 20 percent of SwissProt Eukaryota embeddings TSV file.
    matrix_filepath : str
        Filepath to store the cosine similarity matrix.
    """
    data = pd.read_pickle(embeds_path)
    embeds = np.array(data['embeddings'].tolist(), dtype=np.float32)

    test_data = pd.read_csv(test_accession_path, sep='\t')
    test_accessions = np.array(test_data["accession"].tolist())
    test_data_indices = np.array(test_data["index"].tolist(), dtype=np.int64)

    with ProgressBar(total=test_accessions.shape[0]) as progress:
        similarity_matrix = compute_matrix(embeds, test_data_indices, progress)

    np.savez_compressed(matrix_filepath, similarity_matrix)
    print("Matrix saved")


if __name__ == "__main__":

    # WORD2DOC2VEC
    print("Combination 1")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_1.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_1.npz")

    print("Combination 2")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_2.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_2.npz")

    print("Combination 3")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_3.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_3.npz")

    print("Combination 4")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_4.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_4.npz")

    print("Combination 5")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_5.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_5.npz")

    print("Combination 6")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_6.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_6.npz")

    print("Combination 7")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_7.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_7.npz")

    print("Combination 8")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_8.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_8.npz")

    print("Combination 9")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_9.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_9.npz")

    print("Combination 10")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_10.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_10.npz")

    print("Combination 11")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_11.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_11.npz")

    print("Combination 12")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_12.pkl",
                  "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_word2doc2vec_12.npz")

    # HYBRID
    print("Combination 1")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_1.pkl",
                   "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_1.npz")

    print("Combination 2")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_2.pkl",
                   "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_2.npz")

    print("Combination 3")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_3.pkl",
                   "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_3.npz")

    print("Combination 4")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_4.pkl",
                 "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_4.npz")

    print("Combination 5")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_5.pkl",
                 "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_5.npz")

    print("Combination 6")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_6.pkl",
                 "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_6.npz")

    print("Combination 7")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_7.pkl",
                 "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_7.npz")

    print("Combination 8")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_8.pkl",
                 "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_8.npz")

    print("Combination 9")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_9.pkl",
                 "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_9.npz")

    print("Combination 10")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_10.pkl",
                 "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_10.npz")

    print("Combination 11")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_11.pkl",
                 "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_11.npz")

    print("Combination 12")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_12.pkl",
                 "./data/rev-20220525-UniProtKB-eukaryota-20.tsv", "./data/cosine/cosine_hybrid_12.npz")


