import pandas as pd


def filter_uniref_clusters(input_file: str, output_file: str) -> None:
    """
    Filter entries from Uniref Eukaryota 90% clusters. Removes all cluster only belonging to TrEMBL, only belonging to
    UniParc, only belonging to a combination of TrEMBL and UniParc and writes to a TSV file.
    Parameters
    ----------
    input_file : str
        Filepath to Uniref Eukaryota 90% clusters TSV file.
    output_file : str
        Filepath to filtered Uniref Eukaryota 90% clusters TSV file.
    """
    df = pd.read_csv(input_file, sep='\t')
    print(f"Total number of combinations of entries: {len(df['Types'].unique())}")
    print(df['Types'].unique())
    df.drop(df.index[df['Types'] == 'UniProtKB Unreviewed (TrEMBL)'], inplace=True)
    df.drop(df.index[df['Types'] == 'UniParc'], inplace=True)
    df.drop(df.index[df['Types'] == 'UniProtKB Unreviewed (TrEMBL); UniParc'], inplace=True)
    print(f"Entries after filtering: {len(df)}")
    df.to_csv(output_file, sep='\t')


if __name__ == "__main__":
    filter_uniref_clusters("./data/output/uniref/uniref_ids.tsv", "data/output/uniref/filtered_uniref_cluster_ids.tsv")
