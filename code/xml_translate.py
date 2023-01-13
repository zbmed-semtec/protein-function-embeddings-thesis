import re
import pandas as pd
from typing import Tuple
from bs4 import BeautifulSoup as bs
import xml.etree.ElementTree as ET


def extract_mesh_id(tag: ET.Element) -> str:
    """
    From a matched z:mesh tag, it extracts the correspondent MeSH ID if the
    field "id" is found.
    Parameters
    ----------
    tag : ET.Element
        Object correspondent of a <z:mesh></z:mesh> tag.
    Returns
    -------
    mesh_id: str
        Text with the mesh ID.
    """
    mesh_id_pattern = r"\/MESH\/(.*)"
    if tag.attrib.get("id"):
        mesh_id = "MeSH" + \
                  re.search(mesh_id_pattern, tag.attrib.get("id")).group(1)
    else:
        mesh_id = tag.text.strip()
    return mesh_id


def extract_go_id(tag: ET.Element) -> str:
    """
    From a matched z:go tag, it extracts the correspondent GO ID if the
    field "id" is found.
    Parameters
    ----------
    tag : ET.Element
        Object correspondent of a <z:go></z:go> tag.
    Returns
    -------
    mesh_id: str
        Text with the go ID.
    """
    go_id_pattern = r"\/GO\_(.*)"
    if tag.attrib.get("id"):
        go_id = "GO" + \
                  re.search(go_id_pattern, tag.attrib.get("id")).group(1)
    else:
        go_id = tag.text.strip()
    return go_id


def create_dict(input_filepath: str) -> Tuple[dict, dict]:
    """
    Creates two translation dictionary for MeSH and GO respectively.
    Parameters
    ----------
    input_filepath : str
        File path for the annotated XML file.
    Returns
    -------
    mesh_trans_dict : dict
        MeSH translation dictionary with keys as annotated term and values as its corresponding MeSH ID.
    go_trans_dict : dict
        GO translation dictionary with keys as annotated term and values as its corresponding GO ID.
    """
    namespace = {"z": "https://github.com/zbmed-semtec/protein-function-embeddings-thesis#"}
    mesh_trans_dict = {}
    go_trans_dict = {}
    root = ET.parse(input_filepath).getroot()
    for elem in root.findall('entry', namespace):
        for tagged in elem.findall('function/z:mesh', namespace):
            mesh_id = extract_mesh_id(tagged)
            if not tagged.text in mesh_trans_dict.keys():
                mesh_trans_dict[tagged.text] = mesh_id
        for tagged in elem.findall('function/z:go', namespace):
            go_id = extract_go_id(tagged)
            if not tagged.text in go_trans_dict.keys():
                go_trans_dict[tagged.text] = go_id
    return mesh_trans_dict, go_trans_dict


def translate(input_filepath: str, mesh_trans_dict: dict, go_trans_dict: dict) -> dict:
    """
    Translates the annotated XML file. It looks for the function tag and replaces the annotated terms with its
    corresponding MeSH/GO ID.
    Parameters
    ----------
    input_filepath : str
        File path for annotated XML file.
    mesh_trans_dict : dict
        MeSH translation dictionary with keys as annotated term and values as its corresponding MeSH ID.
    go_trans_dict : dict
        GO translation dictionary with keys as annotated term and values as its corresponding GO ID.
    Returns
    -------
    all_function_text : dict
        Dictionary with keys as Swissprot accession ids and values as function plain text.
    """
    with open(input_filepath, "r") as file:
        # Read each line in the file, readlines() returns a list of lines
        content = file.readlines()
        # Combine the lines in the list into a string
        content = "".join(content)
        bs_content = bs(content, "lxml")
        accession_ids = bs_content.find_all('accession')
        function_text = bs_content.find_all('function')
        all_function_text = {}
        for accession, function in zip(accession_ids, function_text):
            functions = []
            accession_id = accession.text
            for word in function:
                if word.string in mesh_trans_dict:
                    word = mesh_trans_dict[word.string]
                elif word.string in go_trans_dict:
                    word = go_trans_dict[word.string]
                else:
                    word = word.string
                functions.append(word)
            formatted_function = "".join(functions)
            all_function_text[accession_id] = formatted_function
    return all_function_text


def save_text(output_filepath: str, functions_dict: dict) -> None:
    """
    Saves and writes the function plain texts and its corresponding accession ids to a TSV file.
    Parameters
    ----------
    output_filepath : str
        File path for output TSV file.
    functions_dict : dict
        Dictionary with keys as SwissProt accession ids and values as function plain text.
    """
    df = pd.DataFrame(functions_dict.items(), columns=['accession', 'function'])
    df.to_csv(output_filepath, sep='\t')


if __name__ == "__main__":
    mesh_dict, go_dict = create_dict("./data/output/annotations/swissprot_go_mesh_formatted.xml")
    function_plain_text = translate("./data/output/annotations/swissprot_go_mesh_formatted.xml", mesh_dict, go_dict)
    save_text("rev-20220525-UniProtKB-annotated-translated.tsv", function_plain_text)

