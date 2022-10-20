import csv
from gensim.models import Word2Vec
from gensim.parsing.preprocessing import STOPWORDS


def analyze_vocab(model: Word2Vec) -> None:
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


def get_vocab_difference(first_model: Word2Vec, second_model: Word2Vec) -> None:
    """
    Extracts the difference in vocabulary between the two input models.
    Parameters
    ----------
    first_model : Word2Vec
        First WordVec model for comparison.
    second_model : Word2Vec
        Second Word2vec model for comparison.
    """
    first_model = Word2Vec.load(first_model)
    second_model = Word2Vec.load(second_model)
    vocab_one = first_model.wv.key_to_index
    vocab_two = second_model.wv.key_to_index
    difference = list(set(vocab_one).difference(vocab_two))
    cw = csv.writer(open("vocabulary_difference.csv", 'w'))
    cw.writerow(sorted(difference, key=len))
    print("Difference in vocabulary:", len(difference))


if __name__ == "__main__":
    model = Word2Vec.load("./data/output/model/word2doc2vec/cbow/min_count_2/word2vec.model")
    analyze_vocab(model)
    get_vocab_difference("./data/output/model/word2doc2vec/cbow/min_count_2/word2vec.model", "./data/output/model/word2doc2vec/cbow/min_count_3/word2vec.model")