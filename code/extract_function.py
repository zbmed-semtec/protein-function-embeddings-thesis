import re
import csv
import pandas as pd
import xml.etree.ElementTree as ET


def extract_function_comment(input_filepath) -> list[dict[str, str]]:
    """
    Extracts the accession number and the function comment of each entry from the UniProtKB/Swiss-Prot XML file.
    Parameters
    ----------
    input_filepath : str
        Filepath for the UniProtKB/Swiss-Prot XML file.
    Returns
    -------
    all_entries : list[dict[str,str]]
        List of dictionaries of each entry with key as accession number and value as function comment.
    """
    all_entries = []
    root = ET.parse(input_filepath).getroot()
    for elem in root.findall('{http://uniprot.org/uniprot}entry'):
        entry = {}
        accession = elem.find('{http://uniprot.org/uniprot}accession')
        accession_no = accession.text
        entry["Accession"] = accession_no
        comments = elem.findall('{http://uniprot.org/uniprot}comment')
        for comment in comments:
            if comment.attrib['type'] == "function":
                for tag in comment:
                    function = tag.text
                    entry["Function"] = function
            else:
                # Adds an empty string if no function comment is found
                function = ""
        all_entries.append(entry)
    return all_entries


def write_to_tsv(entries, filename):
    """
    Creates and writes each UniPrtKB/Swiss-Prot entry as a row into a tsv file.
    Parameters
    ----------
    entries : list[dict[str,str]
        List of dictionaries of each entry with key as accession number and value as function comment.
    filename : str
        Name for the TSV file.
    """
    tsv_columns = ['Accession', 'Function']
    with open(filename, 'w') as tsv_file:
        writer = csv.DictWriter(tsv_file, fieldnames=tsv_columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(entries)


def preprocess_function(input_tsv_file, output_tsv_file):
    """
    Removes references to the PubMed articles and PMIDs in the function comments and creates a new TSV file.
    Parameters
    ----------
    input_tsv_file : str
        Name of the TSV file containing the accession numbers and the function comments.
    output_tsv_file : str
        Name for the modified TSV file.
    """
    data = pd.read_csv(input_tsv_file, delimiter='\t')
    for index in range(len(data)):
        first_processed_function = re.sub(r" [(]PubMed:[0-9]+(.*?)[)]", '', str(data.loc[index]['Function']))
        processed_function = re.sub(r" [(]Ref.[^)]*[)]", '', str(first_processed_function))
        data.iloc[index, 1] = processed_function
    data.to_csv(output_tsv_file, sep='\t')
    return


if __name__ == "__main__":
    entries = extract_function_comment("data/uniprot_sprot.xml")
    write_to_tsv(entries, "rev-20220525-UniProtKB.tsv")
    preprocess_function("rev-20220525-UniProtKB.tsv", "processed-rev-20220525-UniProtKB.tsv")