import re
import csv
import pandas as pd
import numpy as np
from typing import Any, Iterable
from nltk.tokenize import word_tokenize
from gensim.models import Word2Vec
from gensim.parsing.preprocessing import STOPWORDS


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
        row["Function"] = tokenize(row["Function"])
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


def create_word2vec_model(docs: list, output_file_path: str) -> Word2Vec:
    """
    Generates the Word2Vec model and builds the vocabulary from the function comment texts.
    Parameters
    ----------
    docs : list
        List of functions comments.
    output_file_path : str
        File path for the Word2Vec model.
    Returns
    -------
    model : Word2Vec
        Word2Vec model.
    """
    model = Word2Vec(sentences=docs, vector_size=200, epochs=5, window=5, min_count=4, workers=4)
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
        embeddings_list = []
        for word in functions[index]:
            embeddings_list.append(model.wv[word])
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


def analyze_vocab(model) -> None:
    """
    Outputs the vocabulary statistics and calculates the number of stop words present in the model vocabulary.
    Parameters
    ----------
    model: Word2Vec model
        Gensim word2vec model.
    """
    print("Dimensions of embedding matrix", model.wv.vectors.shape)
    print("Vocabulary size", model.wv.vectors.shape[0])
    vocab = model.wv.key_to_index
    vocab = frozenset(vocab.keys())
    stop_words = STOPWORDS
    print("Number of stop words present in the vocab:", len(vocab.intersection(stop_words)))


def get_vocab_difference(first_model, second_model) -> None:
    first_model = Word2Vec.load(first_model)
    second_model = Word2Vec.load(second_model)
    vocab_one = first_model.wv.key_to_index
    vocab_two = second_model.wv.key_to_index
    difference = list(set(vocab_one).difference(vocab_two))
    cw = csv.writer(open("vocabulary_difference.csv", 'w'))
    cw.writerow(sorted(difference, key=len))
    print("Difference in vocabulary:", len(difference))


if __name__ == "__main__":
    # nltk.download('punkt')
    # accessions, functions = process_from_tsv("./data/output/functions/rev-20220525-UniProtKB.tsv")
    # create_word2vec_model(functions, "./data/output/model/gensim/cbow/min_count_4/word2vec.model")
    # create_document_embeddings(accessions, functions, "./data/output/model/word2vec.model", "./data/output/embeddings")
    # model_sg = Word2Vec.load("./data/output/model/gensim/cbow/min_count_4/word2vec.model")
    # analyze_vocab(model_sg)
    diff = get_vocab_difference("./data/output/model/gensim/sg/min_count_4/word2vec.model", "./data/output/model/gensim/sg/min_count_5/word2vec.model")
    # model_sg_5 = Word2Vec.load("./data/output/model/gensim/sg/min_count_5/word2vec.model")
    # analyze_vocab(model_sg_5)
    # print(model_sg.wv.most_similar("protein"))
    # print(model_cbow.wv.most_similar("protein"))
    # print(model_sg.wv.most_similar("cdh1"))
    # print(model_cbow.wv.most_similar("cdh1"))