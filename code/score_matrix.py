import numpy as np
import pandas as pd


def read_cosine_matrix(matrix_file):
    matrix = np.load(matrix_file)
    cosine_matrix = matrix['arr_0']
    return cosine_matrix


def create_cluster_score_matrix(cluster_file: str, cosine_matrix: np.array, blast_file: str, taxon_file: str, output_file: str) -> None:
    """
    Adds the corresponding cosine similarity score, blast percentage identity score and taxons for each pair of the
    clustered pairs.
    Parameters
    ----------
    cluster_file : str
        TSV Filepath to clustered pairs of proteins.
    cosine_matrix : np.array
        Cosine matrix of dimensions 161354 x 161354 of the best model.
    blast_file : str
        TSV Filepath for blast percentage identity score for all Eukaryota proteins.
    taxon_file : str
        TSV Filepath for Eukaryota fucntions and taxonomy data.
    output_file : str
        Output filepath to save the resulting score matrix matrix.
    """
    df = pd.read_csv(cluster_file, sep='\t')
    ref_indices = df['accession1_index'].to_numpy()
    asd_indices = df['accession2_index'].to_numpy()
    array = np.array(list(zip(ref_indices, asd_indices)), dtype=object)
    # adds cosine score
    cosine_array = np.zeros((len(ref_indices)), dtype=object)
    for idx, i in enumerate(array):
        if i[0] > i[1]:
            cosine_score = cosine_matrix[i[1]][i[0]]
        elif i[0] < i[1]:
            cosine_score = cosine_matrix[i[0]][i[1]]
        cosine_array[idx] = cosine_score
    df['cosine_score'] = cosine_array
    # adds blast percentage identity score
    blast = pd.read_csv(blast_file, sep='\t')
    blast_grouped = blast.groupby(['accession1', 'accession2'])['sequence_identity_score'].max()
    result = pd.merge(df, blast_grouped, on=['accession1', 'accession2'], how='left')
    result['sequence_identity_score'].fillna(
        result.groupby(['accession2', 'accession1'])['sequence_identity_score'].transform('max'), inplace=True)
    # adds taxon for both accessions
    df_euk = pd.read_csv(taxon_file, sep='\t')
    accession_to_taxon = dict(zip(df_euk['accession'], df_euk['taxon']))
    result['taxon_acc_1'] = df['accession1'].map(accession_to_taxon)
    result['taxon_acc_2'] = df['accession2'].map(accession_to_taxon)
    result.to_csv(output_file, sep='\t')
    # # df = pd.read_csv("./data/output/scores/score_hybrid.tsv", sep='\t')
    # # df = df[df['sequence_identity_score'].isnull()]
    # df = pd.read_csv("./add.tsv", sep='\t')
    # print(len(df))
    # df.to_csv("./add2.tsv", sep='\t')
    # result = pd.read_csv("./data/output/scores/clustered_score_matrix_hybrid.tsv", sep='\t')
    # new = pd.read_csv("./values2.tsv", sep='\t')
    # merged = pd.concat([result, new])
    # merged.dropna(subset=['sequence_identity_score'], inplace=True)
    # merged.to_csv("./data/output/scores/clustered_score_matrix_hybrid.tsv", sep='\t')
    # print(len(result.loc[result['cosine_score'] < 0.9]))
    # merged.to_csv("./data/clustered_score_matrix_word2doc2vec.tsv", sep='\t', index=False)


def create_non_clustered_score_matrix(cluster_file: str, cosine_matrix: np.array, output_file: str) -> None:
    """
    Adds the corresponding cosine similarity score for each pair of the
    non-clustered pairs.
    Parameters
    ----------
    cluster_file : str
        TSV Filepath to clustered pairs of proteins.
    cosine_matrix : np.array
        Cosine matrix of dimensions 161354 x 161354 of the best model.
    output_file : str
        Output filepath to save the resulting score matrix matrix.
    """
    df = pd.read_csv(cluster_file, sep='\t')
    ref_indices = df['accession1_index'].to_numpy()
    asd_indices = df['accession2_index'].to_numpy()
    array = np.array(list(zip(ref_indices, asd_indices)), dtype=object)
    cosine_array = np.zeros((len(ref_indices)), dtype=object)
    for idx, i in enumerate(array):
        if i[0] > i[1]:
            cosine_score = cosine_matrix[i[1]][i[0]]
        elif i[0] < i[1]:
            cosine_score = cosine_matrix[i[0]][i[1]]
        cosine_array[idx] = cosine_score
    df['cosine_score'] = cosine_array
    df.to_csv(output_file, sep='\t')


if __name__ == "__main__":
    # Word2doc2Vec
    cosine_matrix = read_cosine_matrix("./data/output/cosine/cosine_word2doc2vev_bestmodel.npz")
    create_cluster_score_matrix("./data/output/uniref/clustered_pairs_index.tsv",
                                cosine_matrix,
                                "./data/output/blast.tsv",
                                "./data/output/functions/rev-20220525-UniProtKB-eukaryota.tsv",
                                "./data/output/scores/clustered_score_matrix_word2doc2vec.tsv")

    create_non_clustered_score_matrix("./data/output/uniref/not_clustered_pairs_index.tsv",
                                      cosine_matrix,
                                      "./data/output/scores/not_clustered_score_matrix_word2doc2vec.tsv")

    # Hybrid-Word2doc2Vec
    cosine_matrix = "./data/output/cosine/cosine_hybrid_bestmodel.npz"
    create_cluster_score_matrix("./data/output/uniref/clustered_pairs_index.tsv",
                                cosine_matrix,
                                "./data/output/blast.tsv",
                                "./data/output/functions/rev-20220525-UniProtKB-eukaryota.tsv",
                                "./data/output/scores/clustered_score_matrix_hybrid.tsv")
    create_non_clustered_score_matrix("./data/output/uniref/not_clustered_pairs_index.tsv",
                                      cosine_matrix,
                                      "./data/output/scores/not_clustered_score_matrix_word2doc2vec.tsv")

