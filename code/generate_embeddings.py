import re
import nltk
import pandas as pd
import numpy as np
import typing
from typing import Any, Iterable
from nltk.tokenize import word_tokenize
from gensim.models import Word2Vec


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


def tokenize(text: str) -> list:
    """
    Converts the text into tokens.
    Parameters
    ----------
    text : str
        Text to be tokenized.
    Returns
    -------
        List of tokens from the text.
    """
    return word_tokenize(text)


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
    while iteration < len(text):
        word = "".join([c for c in str(text[iteration]) if re.match(letters_pattern, c)])
        text[iteration] = word
        iteration += 1
    cleaned_function = []
    for word in text:
        if word != "":
            cleaned_function.append(word)
    return cleaned_function


def preprocess_data(data: pd.DataFrame): # -> tuple[list[Any], list[Any]]:
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
        row["function"] = convert_lowercase(row["function"])
        row["function"] = tokenize(row["function"])
        row["function"] = remove_special_characters(row["function"])

        data.at[index, 'accession'] = row["accession"]
        data.at[index, 'function'] = row["function"]

        accessions.append(row["accession"])
        functions.append(row["function"])
    return accessions, functions


def process_from_tsv(input_file_path: str): # -> tuple[Iterable, Iterable]:
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


def create_word2vec_model(params: dict, docs: list, output_file_path: str): # -> Word2Vec:
    """
    Generates and saves the Word2Vec model and builds the vocabulary from the function comment texts.
    Parameters
    ----------
    params: dict
        Dictionary of input model parameters.
    docs : list
        List of functions comments.
    output_file_path : str
        File path for the Word2Vec model.
    Returns
    -------
    model : Word2Vec
        Word2Vec model.
    """
    params['sentences'] = docs
    model = Word2Vec(**params)
    # model = Word2Vec(sentences=docs, vector_size=200, epochs=5, window=5, min_count=2, workers=4)
    model.build_vocab(docs)
    model.save(output_file_path)
    print("Model saved")
    return model


def create_document_embeddings(accessions: list, functions: list, word2vec_model: Word2Vec, output_dir_path: str) -> None:
    """
    Generates document embeddings from the generated Word2Vec model.
    Parameters
    ----------
    accessions : list
        List of accession numbers.
    functions : list
        List of function comments.
    word2vec_model : Word2Vec
        Generated Word2Vec model.
    output_dir_path: str
        File path for the generated embeddings.
    """
    model = Word2Vec.load(word2vec_model)
    document_embeddings = []

    for index in range(len(accessions)):
        # print(index)
        embeddings_list = []
        for word in functions[index]:
            try:
                embeddings_list.append(model.wv[word])
            # Does not account for words having frequency below min_count.
            except:
                continue
        #  Generate document embeddings from word embeddings
        first = True
        document = []
        for embedding in embeddings_list:
            if first:
                for dimension in embedding:
                    document.append(0.0)
                first = False
            doc_dimension = 0
            for dimension in embedding:
                document[doc_dimension] += dimension
                doc_dimension += 1
        doc_dimension = 0
        for dimension in document:
            # Get the average of each dimension of the embeddings and store it in the document list
            document[doc_dimension] = (dimension / len(embeddings_list))
            doc_dimension += 1
        document_embeddings.append(document)

    for index, pmid in enumerate(accessions):
        np.save(f'{output_dir_path}/{pmid}', document_embeddings[index])

    print("Embeddings Generated")


if __name__ == "__main__":
    # nltk.download('punkt')

    # WORD2DOC2VEC APPROACH
    accessions, functions = process_from_tsv("./data/rev-20220525-UniProtKB-eukaryota.tsv")

    create_document_embeddings(accessions, functions,
                               "./data/model/word2doc2vec/200/cbow/min_count_2/word2vec.model",
                               "./data/embeddings/word2doc2vec/200/cbow/min_count_2/")
    print("Generated model 1")

    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/200/cbow/min_count_3/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/200/cbow/min_count_3/")
    # print("Generated model 2")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/200/cbow/min_count_4/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/200/cbow/min_count_4")
    # print("Generated model 3")

    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/200/sg/min_count_2/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/200/sg/min_count_2/")
    # print("Generated model 4")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/200/sg/min_count_3/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/200/sg/min_count_3/")
    # print("Generated model 5")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/200/sg/min_count_4/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/200/sg/min_count_4/")
    # print("Generated model 6")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/400/cbow/min_count_2/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/400/cbow/min_count_2/")
    # print("Generated model 7")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/400/cbow/min_count_3/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/400/cbow/min_count_3/")
    # print("Generated model 8")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/400/cbow/min_count_4/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/400/cbow/min_count_4/")
    # print("Generated model 9")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/400/sg/min_count_2/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/400/sg/min_count_2/")
    # print("Generated model 10")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/400/sg/min_count_3/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/400/sg/min_count_3/")
    # print("Generated model 11")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/word2doc2vec/400/sg/min_count_4/word2vec.model",
    #                            "./data/embeddings/word2doc2vec/400/sg/min_count_4/")
    # print("Generated model 12")

    # HYBRID APPROACH
    # accessions, functions = process_from_tsv("./data/rev-20220525-UniProtKB-eukaryota-hybrid.tsv")
    #
    # create_document_embeddings(accessions, functions,
    #                           "./data/model/hybridword2doc2vec/200/cbow/min_count_2/word2vec.model",
    #                            "./data/embeddings/hybrid/200/cbow/min_count_2/")
    # print("Generated model 1")
    #
    # create_document_embeddings(accessions, functions,
    #                           "./data/model/hybridword2doc2vec/200/cbow/min_count_3/word2vec.model",
    #                            "./data/embeddings/hybrid/200/cbow/min_count_3/")
    # print("Generated model 2")
    #
    # create_document_embeddings(accessions, functions,
    #                           "./data/model/hybridword2doc2vec/200/cbow/min_count_4/word2vec.model",
    #                            "./data/embeddings/hybrid/200/cbow/min_count_4")
    # print("Generated model 3")
    #
    # create_document_embeddings(accessions, functions,
    #                           "./data/model/hybridword2doc2vec/200/sg/min_count_2/word2vec.model",
    #                            "./data/embeddings/hybrid/200/sg/min_count_2/")
    # print("Generated model 4")
    #
    # create_document_embeddings(accessions, functions,
    #                           "./data/model/hybridword2doc2vec/200/sg/min_count_3/word2vec.model",
    #                            "./data/embeddings/hybrid/200/sg/min_count_3/")
    # print("Generated model 5")
    #
    # create_document_embeddings(accessions, functions,
    #                           "./data/model/hybridword2doc2vec/200/sg/min_count_4/word2vec.model",
    #                            "./data/embeddings/hybrid/200/sg/min_count_4/")
    # print("Generated model 6")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/hybridword2doc2vec/400/cbow/min_count_2/word2vec.model",
    #                            "./data/embeddings/hybrid/400/cbow/min_count_2/")
    # print("Generated model 7")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/hybridword2doc2vec/400/cbow/min_count_3/word2vec.model",
    #                            "./data/embeddings/hybrid/400/cbow/min_count_3/")
    # print("Generated model 8")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/hybridword2doc2vec/400/cbow/min_count_4/word2vec.model",
    #                            "./data/embeddings/hybrid/400/cbow/min_count_4/")
    # print("Generated model 9")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/hybridword2doc2vec/400/sg/min_count_2/word2vec.model",
    #                            "./data/embeddings/hybrid/400/sg/min_count_2/")
    # print("Generated model 10")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/hybridword2doc2vec/400/sg/min_count_3/word2vec.model",
    #                            "./data/embeddings/hybrid/400/sg/min_count_3/")
    # print("Generated model 11")
    #
    # create_document_embeddings(accessions, functions,
    #                            "./data/model/hybridword2doc2vec/400/sg/min_count_4/word2vec.model",
    #                            "./data/embeddings/hybrid/400/sg/min_count_4/")
    # print("Generated model 12")
