import pandas as pd


def group_by_taxon(file_path: str) -> None:
    """
    Analyzes different levels of kingdoms and returns count of entries based on each kingdom.
    Parameters
    ----------
    file_path : str
        File path of the TSV containing metadata of SwissProt entries.
    """
    data = pd.read_csv(file_path, sep='\t')
    kingdom_level = []
    eukaryota_count = 0
    archaea_count = 0
    virus_count = 0
    bacteria_count = 0
    for idx, row in data.iterrows():
        kingdom_level_value = row['taxonomy_lineage'].split(',')[0].strip('[').strip(']')
        kingdom_level.append(kingdom_level_value)
        if kingdom_level_value == "'Eukaryota'":
            eukaryota_count += 1
        elif kingdom_level_value == "'Archaea'":
            archaea_count += 1
        elif kingdom_level_value == "'Viruses'":
            virus_count += 1
        else:
            bacteria_count += 1
    taxon_level = data['taxon'].nunique()
    print("Kingdom level classification", len(set(kingdom_level)))
    print("Kingdoms:", set(kingdom_level))
    print("Eukaryota:", eukaryota_count)
    print("Archaea:", archaea_count)
    print("Viruses:", virus_count)
    print("Bacteria:", bacteria_count)
    print("Taxon level classification", taxon_level)


if __name__ == "__main__":
    group_by_taxon("./data/output/functions/rev-20220525-UniProtKB.tsv")