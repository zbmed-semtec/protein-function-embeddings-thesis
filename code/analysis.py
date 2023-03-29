import pandas as pd


def analyse_non_clusters(uniref_id_mapping_file: str, non_clustered_file: str, taxon_file: str) -> None:
    """
    Analyses and compares the non clustered pairs with the latest release 2023_01 of UniProt to check if any of the
    pairs have already been clustered. Adds this information as a new column "clustered", sets it to True is clustered.
    Additionally, adds Taxon name of both the accessions of each pair as two additonal columns and creates a file
    containing of all pairs that are to date still not clustered (true non clusters).
    Parameters
    ----------
    uniref_id_mapping_file : str
        TSV file consisting of information about each Eukaryota entry and if they belong to any cluster
        as per release 2023_01.
    non_clustered_file : str
        TSV filepath to score matrix of Non clustered file.
    taxon_file : str
        TSV Filepath for Eukaryota fucntions and taxonomy data.
    """
    df_id_mapping = pd.read_csv(uniref_id_mapping_file, sep='\t', usecols=['From', 'Cluster ID'], index_col=0, squeeze=True)
    dictionary = {index: value for index, value in df_id_mapping.items()}

    df1 = pd.read_csv(non_clustered_file, sep='\t')
    # adds True if the pair has been clustered as per the latest release, else False
    for idx, row in df1.iterrows():
        if dictionary[row['accession1']] == dictionary[row['accession2']]:
            df1.at[idx, 'clustered'] = True
        else:
            df1.at[idx, 'clustered'] = False

    # Saves along with an additional 'clustered' column
    df1.to_csv("./data/output/scores/not_clustered_score_matrix_word2doc2vec_with_new_clusters.tsv", sep='\t')

    # adds taxon for both accessions
    df_euk = pd.read_csv(taxon_file, sep='\t')
    accession_to_taxon = dict(zip(df_euk['accession'], df_euk['taxon']))
    df1['taxon_acc_1'] = df1['accession1'].map(accession_to_taxon)
    df1['taxon_acc_2'] = df1['accession2'].map(accession_to_taxon)
    df1.to_csv("./data/output/scores/not_clustered_score_matrix_word2doc2vec_with_taxon.tsv", sep='\t')

    # Saves only those that have not been clustered
    df1 = df1.loc[df1['clustered'] == True]
    df1.to_csv("./data/output/scores/true_non_clusters.tsv", sep='\t')


def analyse_cases(true_non_cluster_file: str):
    """
    Analysis the true non clusters and groups them into two cases and creates two separate TSV files.
    Case 1: high sequence identity and high cosine similarity
    Case 2:  high sequence identity and low cosine similarity.
    Parameters
    ----------
    true_non_cluster_file : str
        TSV Filepath containing of all true non clusters.
    """

    df = pd.read_csv(true_non_cluster_file, sep='\t')
    df1 = df.loc[(df['sequence_identity_score'] >= 0.90) & (df['cosine_score'] >= 0.90)]
    df1.to_csv("./data/output/score/case1.tsv", sep='\t')
    df2 = df.loc[(df['sequence_identity_score'] >= 0.90) & (df['cosine_score'] < 0.90)]
    df2.to_csv("./data/output/score/case2.tsv", sep='\t')


if __name__ == "__main__":
    analyse_non_clusters("./data/output/scores/uniprot-compressed_true_download_true_fields_id_2Cname_2Ctypes_2Ccou-2023.03.09-14.53.27.42.tsv",
                         "./data/output/scores/not_clustered_score_matrix_word2doc2vec.tsv", "./data/output/functions/rev-20220525-UniProtKB-eukaryota.tsv")
    analyse_cases("./data/output/scores/true_non_clusters.tsv")
