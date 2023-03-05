import csv
import itertools
import pandas as pd


def read_cluster_index(cluster_index_file: str) -> pd.DataFrame:
    """
    Reads the cluster index file containing Eukaryota accessions and its corresponding cluster index.
    Parameters
    ----------
    cluster_index_file : str
        Filepath to the cluster index file.

    Returns
    -------
    df : pd.Dataframe
        Pandas dataframe of the cluster index file. Consists of four columns: accession_index, accession,
        tp/fp and cluster_index.
    """
    df = pd.read_csv(cluster_index_file, sep='\t', names=['accession_index', 'accession', 'tp/fp', 'cluster_index'])
    to_array = lambda x: eval(x)
    df['cluster_index'] = df['cluster_index'].apply(to_array)
    return df


def group_clusters(cluster_indices: pd.DataFrame, cluster_filename: str, non_cluster_filename: str) -> None:
    """
    Groups Eukaryota accessions into two dataframes: those belonging to a cluster and those not belonging to any cluster.
    Saves both the dataframes as TSV files.
    Parameters
    ----------
    cluster_indices : pd.Dataframe
        Pandas dataframe of the cluster index file.
    cluster_filename : str
        Filepath to save accessions belonging to a cluster.
    non_cluster_filename : str
        Filepath to save accessions not belonging to a cluster.
    """
    # all accessions belonging to a cluster
    clusters = cluster_indices[cluster_indices['tp/fp'] == 1]
    clusters.to_csv(cluster_filename, sep='\t', index=False)

    # all accessions not belonging to a cluster
    not_clusters = cluster_indices[cluster_indices['tp/fp'] == 0]
    not_clusters.to_csv(non_cluster_filename, sep='\t', index=False)


def generate_cluster_pairs(clusters: str, output_filename: str) -> None:
    """
    Generates accession pairs between all accessions belonging to the same cluster and saves them as a TSV file.
    Parameters
    ----------
    clusters : str
        Filepath of accessions belonging to a cluster.
    output_filename : str
        Filepath to save the pairs belonging to a cluster.
    """
    df = pd.read_csv(clusters, sep='\t')
    to_array = lambda x: eval(x)
    df['cluster_index'] = df['cluster_index'].apply(to_array)
    df = df.sort_values(by='cluster_index')
    df['cluster_index'] = df['cluster_index'].apply(lambda x: x[0])
    main_dict = dict()
    for idx, row in df.iterrows():
        accession = row['accession']
        index = row['cluster_index']
        frequency = df[df['cluster_index'] == index].shape[0]
        if frequency >= 2:
            entries = df['accession'].loc[df['cluster_index'] == index].to_list()
            main_dict[index] = entries
    print("Processed counts")
    long_list = []
    for key, value in main_dict.items():
        combinations = itertools.combinations(value, 2)
        for combination in combinations:
            long_list.append([combination, key])
    print(len(long_list))
    print("Processed combinations")

    with open(output_filename, 'w') as op:
        writer = csv.writer(op, delimiter='\t')
        for row in long_list:
            row = list(row[0]) + [row[1]]
            writer.writerow(row)


def generate_non_cluster_pairs(non_clusters: str, blast_file: str, output_filename: str) -> None:
    """
    Generates accession pairs between all accessions not belonging to any cluster. Reads the non clusters and blast
    result file, extracts the BLAST results for each of the accession, creates a pair and stores the BLAST percentage
    identity score for the pair.
    Parameters
    ----------
    non_clusters : str
        Filepath of accessions not belonging to a cluster.
    blast_file : str
        Filepath of blast results for all Eukaryota accessions.
    output_filename : str
        Filepath to save the pairs not belonging to a cluster.
    """
    not_clusters = pd.read_csv(non_clusters, sep='\t')
    accessions = not_clusters['accession'].to_list()
    blast = pd.read_csv(blast_file, sep='\t')
    entries = blast[blast.isin(accessions).any(axis=1)]
    entries = entries[entries['sequence_identity_score'] >= 90.0]
    entries = entries[entries['accession1'] != entries['accession2']]
    entries.to_csv(output_filename, sep='\t', index=None)


def add_index(file: str, output_filename: str, clustered: bool = True) -> None:
    """
    Adds the index value for the both the accessions of a pair and saves it as a TSV file.
    Parameters
    ----------
    file : str
        Filepath to either pairs belonging or not belonging to a cluster.
    output_filename  :str
        Filepath to save the file along with the indices.
    clustered : bool
        Checks if the file belongs to clustered pairs or non-clustered pairs. Default value is True.
    """
    if clustered:
        df = pd.read_csv(file, sep='\t', names=['accession1', 'accession2', 'cluster_index'])
    else:
        df = pd.read_csv(file, sep='\t')
    data = pd.read_pickle("./data/output/embeddings_pickle/word2doc2vec/word2doc2vec_embs_euk_4.pkl")
    accession_to_index = dict(zip(data['accessions'], data.index))

    for index, row in df.iterrows():
        try:
            value1 = int(accession_to_index.get(row['accession1']))
            df.at[index, 'accession1_index'] = value1
        except:
            df = df.drop(index)
        try:
            value2 = int(accession_to_index.get(row['accession2']))
            df.at[index, 'accession2_index'] = value2
        except:
            df = df.drop(index)

    df = df.astype({'accession1_index': int})
    df = df.astype({'accession2_index': int})
    df.to_csv(output_filename, sep='\t')


if __name__ == "__main__":
    cluster_indices = read_cluster_index("./data/output/uniref/cluster_index.tsv")
    group_clusters(cluster_indices, "./data/output/uniref/clusters.tsv", "./data/output/uniref/not_clusters.tsv")
    generate_cluster_pairs("./data/output/uniref/clusters.tsv", "./data/output/uniref/clustered_pairs.tsv")
    add_index("./data/output/uniref/clustered_pairs.tsv", "./data/output/uniref/clustered_pairs_index.tsv")
    generate_non_cluster_pairs("./data/output/uniref/not_clusters.tsv", "./data/output/blast.tsv",
                               "./data/output/uniref/not_clustered_pairs.tsv")
    add_index("./data/output/uniref/not_clustered_pairs.tsv",
              "./data/output/uniref/not_clustered_pairs_index.tsv", clustered=False)
