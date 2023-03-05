import re
import csv
import itertools
import pandas as pd
import xml.etree.ElementTree as ET


def extract_function_comment(input_filepath) -> list[dict[str, str]]:
    """
    Extracts the accession number, function comment, evidence tags, taxon, taxonomy lineage and taxon id of each entry
    from the UniProtKB/Swiss-Prot XML file.
    Parameters
    ----------
    input_filepath : str
        Filepath for the UniProtKB/Swiss-Prot XML file.
    Returns
    -------
    all_entries : list[dict[str,str]]
        List of dictionaries of each entry with key as accession number and value as function comment, evidence tags,
        and taxonomy metadata.
    """
    all_entries = []
    root = ET.parse(input_filepath).getroot()
    for elem in root.findall('{http://uniprot.org/uniprot}entry'):
        entry = {}
        accession = elem.find('{http://uniprot.org/uniprot}accession')
        accession_no = accession.text
        entry["accession"] = accession_no
        comments = elem.findall('{http://uniprot.org/uniprot}comment')
        for comment in comments:
            if comment.attrib['type'] == "function":
                for tag in comment:
                    function = tag.text
                    entry["function"] = function
                    if 'evidence' in tag.attrib:
                        # All evidence tags in the function comment.
                        evidence_tags = tag.attrib['evidence'].split(" ")
                        all_eco_tags = []
                        # All evidence tags in the XML file of the protein.
                        evidence = elem.findall('{http://uniprot.org/uniprot}evidence')
                        for value in evidence_tags:
                            eco_tag = [tag.attrib['type'] for tag in evidence if tag.attrib['key'] == value]
                            all_eco_tags.append(eco_tag[0])
                        entry["evidence_tags"] = ", ".join(all_eco_tags)
            else:
                # Adds an empty string if no function comment is found.
                function = ""
        taxonomy = elem.findall('{http://uniprot.org/uniprot}organism')
        for data in taxonomy:
            for subtag in data:
                if subtag.tag == "{http://uniprot.org/uniprot}dbReference":
                    taxon_id = subtag.attrib["type"] + " " + subtag.attrib["id"]
                if subtag.tag == "{http://uniprot.org/uniprot}lineage":
                    taxonomy_lineage = []
                    for lineage in subtag:
                        taxonomy_lineage.append(lineage.text)
                    entry["taxon"] = taxonomy_lineage[-1]
                    entry["taxonomy_lineage"] = taxonomy_lineage
                    entry["taxon_identifier"] = taxon_id
        all_entries.append(entry)
    return all_entries


def write_to_tsv(entries, filename):
    """
    Creates and writes each UniProtKB/Swiss-Prot entry as a row into a tsv file.
    Parameters
    ----------
    entries : list[dict[str,str]
        List of dictionaries of each entry with key as accession number and value as function comment.
    filename : str
        Name for the TSV file.
    """
    tsv_columns = ['accession', 'function', 'evidence_tags', 'taxon', 'taxonomy_lineage', 'taxon_identifier']
    with open(filename, 'w') as tsv_file:
        writer = csv.DictWriter(tsv_file, fieldnames=tsv_columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(entries)


def preprocess_function(input_tsv_file, output_tsv_file):
    """
    Removes references to the PubMed articles, PMIDs and evidence tags in the function comments and creates a new TSV file.
    Parameters
    ----------
    input_tsv_file : str
        Name of the TSV file containing the accession numbers and the function comments.
    output_tsv_file : str
        Name for the modified TSV file.
    """
    data = pd.read_csv(input_tsv_file, delimiter='\t', dtype='str')
    # Adds a column for all reference texts.
    data.insert(3, "evidence_text", "")
    # Drops all entries without a function.
    data = data.dropna(subset=['Function']).reset_index(drop=True)
    for index in range(len(data)):
        pubmed_refs = re.findall(r"[(]PubMed:[0-9]+[)]", str(data.loc[index]["Function"]))  # PubMed Articles
        refs = re.findall(r"[(]Ref.[^)]*[)]", str(data.loc[index]["Function"]))  # References to PMIDs
        sim_tags = re.findall(r"[(]By similarity[)]", str(data.loc[index]["Function"]))  # Evidence tags
        evidence_text = list(itertools.chain(pubmed_refs, refs, sim_tags))
        first_processed_function = re.sub(r" [(]PubMed:[0-9]+(.*?)[)]", '', str(data.loc[index]['Function']))  # PubMed Articles
        second_processed_function = re.sub(r" [(]Ref.[^)]*[)]", '', str(first_processed_function))  # References to PMIDs
        processed_function = re.sub(r" [(]By similarity[)]", '', str(second_processed_function))  # Evidence tags
        data.iloc[index, 1] = processed_function
        data.iloc[index, 3] = ", ".join(evidence_text)
    data.to_csv(output_tsv_file, sep='\t')
    return


if __name__ == "__main__":
    entries = extract_function_comment("data/uniprotkb/reviewed/uniprot_sprot.xml")
    write_to_tsv(entries, "rev-20220525-UniProtKB.tsv")
    preprocess_function("data/output/functions/rev-20220525-UniProtKB.tsv", "data/processed-rev-20220525-UniProtKB.tsv")