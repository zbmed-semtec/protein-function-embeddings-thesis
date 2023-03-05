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
def compute_matrix(embeds: np.array, progress_proxy: tqdm) -> np.array:
    """
    Creates an empty numpy matrix with rows and columns based on the dimensions of embeddings. Iterates through each
    cell of the matrix to calculate the cosine similarity between the pair of embeddings. Stores the cosine similarity
    score and returns the filled matrix.
    Parameters
    ----------
    embeds : np.array
        All SwissProt Eukaryota embeddings in a numpy array.
    progress_proxy : tqdm
        tqdm progress bar proxy object to show progress.
    Returns
    -------
    matrix : Cosine similarity matrix.
    """
    matrix = np.zeros((embeds.shape[0], embeds.shape[0]), dtype=np.float32)
    for i in range(embeds.shape[0]):
        embed1 = embeds[i]
        for j in range(embeds.shape[0]):
            embed2 = embeds[j]
            cs = round(cosine_similarity_numba(embed1, embed2), 4)
            if j > i:
                matrix[i, j] = cs
        progress_proxy.update(1)
    return matrix


def create_matrix(embeds_path: str, matrix_filepath: str) -> None:
    """
    Wrapper function to read the SwissProt Eukaryota embeddings pickle file, compute cosine similarity and save the
    cosine similarity matrix as a .NPZ file.
    Parameters
    ----------
    embeds_path : str
        Filepath to SwissProt Eukaryota embeddings pickle file.
    matrix_filepath : str
        Filepath to store the cosine similarity matrix.
    """
    data = pd.read_pickle(embeds_path)
    embeds = np.array(data['embeddings'].tolist(), dtype=np.float32)

    with ProgressBar(total=embeds.shape[0]) as progress:
        similarity_matrix = compute_matrix(embeds, progress)

    np.savez_compressed(matrix_filepath, similarity_matrix)
    print("Matrix saved")


if __name__ == "__main__":
    print("Best model word2doc2vec")
    create_matrix("./data/emb_pickle/word2doc2vec_embs_euk_4.pkl",
                  "./data/best_model/cosine_word2doc2vec_bestmodel.npz")

    print("Best model hybrid-word2doc2vec")
    create_matrix("./data/emb_pickle/hybrid_embs_euk_4.pkl",
                  "./data/best_model/cosine_hybrid_bestmodel.npz")

