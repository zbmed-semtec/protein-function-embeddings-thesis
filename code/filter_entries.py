import pandas as pd


def extract_superkingdom(input_filename: str, output_filename: str) -> None:
    """
    Reads the SwissProt all entries TSV file, extracts the superkingdom and adds it as a separate column. Writes this
    to a new TSV file.
    Parameters
    ----------
    input_filename : str
        Filepath to the SwissProt all entries TSV file.
    output_filename : str
        Filepath to the SwissProt all entries with superkingdom TSV file.
    """
    data = pd.read_csv(input_filename, sep='\t', index_col=0)
    data.insert(5, "kingdom", "")
    for idx, row in data.iterrows():
        kingdom_level_value = row['taxonomy_lineage'].split(',')[0].strip('[').strip(']').strip("'")
        data.iloc[idx, 5] = kingdom_level_value
    data.to_csv(output_filename, sep='\t')


def filter_entries(input: str, output: str) -> None:
    """
    Reads the SwissProt all entries TSV file with superkingdom, filters out all Eukaryota entries and writes it to a new
    TSV file.
    Parameters
    ----------
    input : str
        Filepath to the SwissProt all entries with superkingdom TSV file.
    output : str
        Filepath to the SwissProt Eukaryota entries TSV file.
    """
    df = pd.read_csv(input, sep='\t', index_col=0)
    new_df = df.loc[df['kingdom'] == "Eukaryota"]
    new_df = new_df.reset_index(drop=True)
    new_df.to_csv(output, sep='\t')


def filter_test_eukaryota(input: str, percent: float, output: str) -> None:
    """
    Reads the SwissProt Eukaryota TSV file and extracts the given input percent of entries from it and writes it to a
    new TSV file.
    Parameters
    ----------
    input : str
        Filepath to the SwissProt Eukaryota entries TSV file.
    percent : float
        Percentage of entries to filer.
    output : str
        Filepath to the SwissProt Eukaryota 20 percent entries TSV file.
    """
    df = pd.read_csv(input, sep='\t', index_col=0)
    df.sort_values(by="accession", ignore_index=True, inplace=True)
    test_sample = df.sample(frac=percent)
    test_sample.sort_values(by="accession", inplace=True)
    test_sample = test_sample.reset_index()
    test_sample.to_csv(output, sep='\t')


def filter_hybrid_entries(input_eukaryota_entries: str, input_hybrid_eukaryota: str, output_hybrid_eukaryota: str) -> None:
    """
    Reads the SwissProt Eukaryota TSV file and extracts those accession entries from the annotated SwissProt all entries
    TSV file writes it to a new TSV file.
    Parameters
    ----------
    input_eukaryota_entries : str
        Filepath to the SwissProt Eukaryota entries TSV file.
    input_hybrid_eukaryota : str
        Filepath to the SwissProt Eukaryota annotated entries TSV file.
    output_hybrid_eukaryota : str
        Filepath to the SwissProt Eukaryota hybrid entries TSV file.
    """
    df = pd.read_csv(input_eukaryota_entries, sep='\t', index_col=0)
    df_hybrid = pd.read_csv(input_hybrid_eukaryota, sep='\t', index_col=0)
    merged_df = df_hybrid.merge(df, on='accession')
    filtered_df = merged_df[merged_df['accession'].notnull()]
    filtered_df = filtered_df[['accession', 'function_x']]
    filtered_df = filtered_df.rename(columns={'accession': 'accession', 'function_x': 'function'})
    filtered_df.to_csv(output_hybrid_eukaryota, sep='\t')


if __name__ == "__main__":
    extract_superkingdom("./data/output/functions/rev-20220525-UniProtKB.tsv",
                         "./data/output/functions/rev-20220525-UniProtKB-kingdoms.tsv")
    filter_entries("./data/output/functions/rev-20220525-UniProtKB-kingdoms.tsv",
                   "./data/output/functions/rev-20220525-UniProtKB-eukaryota.tsv")
    filter_test_eukaryota("./data/output/functions/rev-20220525-UniProtKB-eukaryota.tsv", 0.2,
                          "./data/output/functions/rev-20220525-UniProtKB-eukaryota-20.tsv")
    filter_hybrid_entries("./data/output/functions/rev-20220525-UniProtKB-eukaryota.tsv",
                          "./data/output/functions/rev-20220525-UniProtKB-annotated-translated.tsv",
                          "./data/output/functions/rev-20220525-UniProtKB-eukaryota-hybrid.tsv")
