import pandas as pd
import matplotlib.pyplot as plt


def plot_sequence_embedding(input_filepath: str, output_filepath: str):
    """
    Plots a scatter plot between sequence (blast percentage identity score) vs embedding similarity (cosine similarity)
    for clustered pairs.
    Parameters
    ----------
    input_filepath : str
        TSV filepath for the clustered pairs of proteins.
    output_filepath : str
        Filepath to save the plot.
    """
    df = pd.read_csv(input_filepath, sep='\t')
    x = df['cosine_score']
    y = df['sequence_identity_score']
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y, s=5)
    ax.plot([min(x), max(x)], [min(y), max(y)], color='red')
    ax.set_xlabel('Embedding similarity (cosine similarity)')
    ax.set_ylabel('Sequence similarity (percentage identity)')
    ax.set_title('Sequence vs. Embedding similarity')
    plt.savefig(output_filepath)


if __name__ == "__main__":
    plot_sequence_embedding("./data/output/scores/clustered_score_matrix_word2doc2vec.tsv",
                            './data/output/plots/scatter_plot_not_clustered_word2odc2vec.png')

    plot_sequence_embedding("./data/output/scores/clustered_score_matrix_hybrid.tsv",
                            './data/output/plots/scatter_plot_not_clustered_hybrid.png')