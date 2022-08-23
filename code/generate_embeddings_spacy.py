import os
import re
import pandas as pd
import numpy as np
import spacy
from typing import Any, Iterable


def convert_lowercase(text: str) -> str:
    """
    Converts the text into lowercase
    Parameters
    ----------
    text : str
        Text to be converted.
    Returns
    -------
        Text converted into lowercase.
    """
    return str(text).lower()


def remove_special_characters(text: str):
    """
    Removes all special characters from the text and only keeps all alphabets, numbers and hyphen.
    Parameters
    ----------
    text : str
        Text to be processed.
    Returns
    -------
        Processed text without special characters.
    """
    letters_pattern = r'.*[a-zA-Z\d\-].*'  # Includes all letters which are numbers, letters or a hyphen.
    iteration = 0
    text = text.split(" ")
    while iteration < len(text):
        word = "".join([c for c in str(text[iteration]) if re.match(letters_pattern, c)])
        text[iteration] = word
        iteration += 1
    cleaned_function = []
    for word in text:
        if word != "":
            cleaned_function.append(word)
    return " ".join(cleaned_function)


def preprocess_data(data: pd.DataFrame) -> tuple[list[Any], list[Any]]:
    """
    Wrapper function for the pre-processing of the text.
    Parameters
    ----------
    data : pd.DataFrame
        Dataframe containing accession numbers and function comments.
    Returns
    -------
    acccessions : list
        List of all accession numbers.
    functions: list
        List of processed function comments.
    """
    accessions = []
    functions = []
    for index, row in data.iterrows():
        row["Function"] = convert_lowercase(row["Function"])
        row["Function"] = remove_special_characters(row["Function"])
        data.at[index, 'Accession'] = row["Accession"]
        data.at[index, 'Function'] = row["Function"]

        accessions.append(row["Accession"])
        functions.append(row["Function"])
    return accessions, functions


def process_from_tsv(input_file_path: str) -> tuple[Iterable, Iterable]:
    """
    Loads and pre-processes the data from the input TSV file.
    Parameters
    ----------
    input_file_path : str
        Filepath for the input TSV file.
    """
    data = pd.read_csv(input_file_path, sep='\t')
    accessions, functions = preprocess_data(data)
    return accessions, functions


def load_model():
    """Loads and returns the downloaded spacy en_core_eng_md model."""
    nlp = spacy.load('en_core_web_md')
    return nlp


def create_embeddings(accessions: list, functions: list, nlp: spacy.model, output_dir_path: str):
    """
    Generates document embeddings from the generated SpaCy model.
    Parameters
    ----------
    accessions : list
        List of accession numbers.
    functions : list
        List of function comments.
    nlp :
        Generated SpaCy model.
    output_dir_path : str
        File path for the generated embeddings.
    """
    document_embeddings = []
    for index in range(len(accessions)):
        # Generate word embeddings
        embeddings_list = []
        doc = nlp(functions[index])
        for word_index in range(len(doc)):
            embeddings_list.append(doc[word_index].vector)
        # Generate document embeddings
        document_embeddings.append(doc.vector)

    for index, pmid in enumerate(accessions):
        np.save(f'{output_dir_path}/{pmid}', document_embeddings[index])


if __name__ == "__main__":
    # os.sys("python3 -m spacy download en_core_web_md")
    accessions, functions = process_from_tsv("./data/output/functions/rev-20220525-UniProtKB.tsv")
    nlp = load_model()
    create_embeddings(accessions, functions, nlp, "./data/output/embeddings")