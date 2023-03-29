import pandas as pd
import numpy as np
from Bio import SeqIO


def extract_eukaryota_accessions(eukaryota_file: str) -> np.array:
    """
    Extracts all accessions from the Eukaryota file and returns them as a numpy array.
    Parameters
    ----------
    eukaryota_file : str
        Filepath to Eukaryota accessions.
    Returns
    -------
    array : np.array
        Numpy array of all accessions.
    """
    file = pd.read_csv(eukaryota_file, sep='\t')
    array = file['accession'].to_numpy()
    return array


def parse_fasta(input_fasta_file: str, eukaryota_accessions: np.array) -> None:
    """
    Extracts all sequences of Eukaryota entries from the SwissProt fasta file and writes them into a fasta file.
    Parameters
    ----------
    input_fasta_file : str
        Filepath to SwissProt fasta file.
    eukaryota_accessions : np.array
        Numpy array of all Eukaryota accessions.
    """
    sequences = [i for i in SeqIO.parse(input_fasta_file, 'fasta')]
    with open("swissprot-eukaryota.fasta", "a") as output:
        for entry in sequences:
            seq_name = entry.name
            accession = seq_name.split("|")[1]
            if accession in eukaryota_accessions:
                SeqIO.write(entry, output, "fasta")


def batch_iterator(iterator: int, batch_size: int) -> list:
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.Align.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    Parameters
    ----------
    iterator : int
        Iterator object
    batch_size : int
        Length of batch size.
    Returns
    ----------
    batch : list
        Generator object.
    """
    for i in range(0, len(iterator), batch_size):
        batch = iterator[i:i + batch_size]
        yield batch


def write_fasta(fast_file: str) -> None:
    """
    Splits the fasta file into batches of 10,000 sequences per file. Creates 17 fasta files in total.
    Parameters
    fasta_file : str
        Filepath to the fasta file containing all Eukaryota sequences.
    """
    record_iter = [i for i in SeqIO.parse(open(fast_file), "fasta")]
    for i, batch in enumerate(batch_iterator(record_iter, 10000)):
        filename = "eukaryota_group_%i.fasta" % (i + 1)
        with open(filename, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, filename))


if __name__ == "__main__":
    eukaryota_accessions = extract_eukaryota_accessions("./data/output/functions/rev-20220525-UniProtKB-eukaryota.tsv")
    eukaryota_fasta = parse_fasta("./data/uniprot/swissprot/uniprot_sprot-only2022_02/uniprot_sprot.fasta", eukaryota_accessions)
    write_fasta("./data/output/functions/swissprot-eukaryota.fasta")


