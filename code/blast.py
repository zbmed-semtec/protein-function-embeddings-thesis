import os
import csv


def merge_blast_files(input_root_directory: str) -> None:
    """
    Merges and uncompressed all batch-wise blast result files into one putput file.
    Parameters
    ----------
    input_root_directory : str
        Input directory containing all blast results files.
    """
    os.chdir("./blast")
    for subdir in os.listdir(input_root_directory):
        if subdir.startswith("Group"):
            os.system(f"cat *.out.gz > ./blast/files/{subdir}.out.gz")
    os.chdir("./blast/files/")
    os.system("cat *.out.gz > blast.out.gz")
    os.system("gzip -d blast.out.gz")


def parse_blast_results(blast_file: str) -> list:
    """
    Parses blast result file and generates pairs of accessions along with the blast identity score.
    Returns these values as as list.
    Parameters
    ----------
    blast_file : str
        Filepath of blast result file.
    Returns
    -------
    blast_results : list
        Nested list containing of accession1, accession2 and sequence identity score.
    """
    blast_results = []
    with open(blast_file, 'r') as file:
        content = file.readlines()
        for line in content:
            if line.startswith('sp|'):
                line_content = line.split("\t")
                key = line_content[0].split("|")[1]
                value = line_content[1].split(".")[0]
                score = line_content[2]
                blast_results.append([key, value, score])

    return blast_results


def blast_to_tsv(blast_results: list, output_blast_file: str) -> None:
    """
    Writes the blast pair results to a TSV file.
    Parameters
    ----------
    blast_results : list
        Nested list containing of accession1, accession2 and sequence identity score.
    output_blast_file : str
        Filepath to save the blast TSV file.
    """
    with open(output_blast_file, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['accession1', 'accession2', 'sequence_identity_score'])
        writer.writerows(blast_results)


if __name__ == "__main__":
    merge_blast_files("./blast")
    results = parse_blast_results("./blast.out")
    blast_to_tsv(results, "./blast.tsv")
