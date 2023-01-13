import gzip
import requests
import pandas as pd


def request_api(url: str, filename: str) -> None:
    """
    Sends a REST API request to UniProt's API to extract Uniref Eukaryota clusters belonging to 90% identity with a cluster size
    greater than or equal to 2 using the streaming endpoint.
    Parameters
    ----------
    url : str
        API URL.
    filename : str
        Filepath to store the TSV data.
    """

    response = requests.get(url)
    print(response.status_code)
    content = gzip.decompress(response.content).decode('utf-8')
    with open(filename, 'w') as write_file:
        write_file.write(content)


def merge_outputs(first_file: str, second_file: str, output_file: str) -> None:
    """
    Merges both the TSV files into one TSV containing all Uniref 90% Eukaryota >= 2 member cluster IDs.
    Parameters
    ----------
    first_file : str
        First file containing about 5,436,578 entries.
    second_file : str
        Second file containing about 5,839,442 entries.
    output_file : str
        Filepath to save the merged TSV file containing about 11,276,020 entries.
    """
    df1 = pd.read_csv(first_file, sep='\t')
    df2 = pd.read_csv(second_file, sep='\t')
    pd.concat([df1, df2]).to_csv(output_file, sep='\t')
    df = pd.read_csv(output_file, sep='\t')
    print(len(df))



if __name__=="__main__":
    url1 = "https://rest.uniprot.org/uniref/stream?compressed=true&fields=id%2Cname%2Ctypes%2Ccount&format=tsv&query=%28%28taxonomy_id%3A2759%29%20AND%20%28count%3A%5B2%20TO%202%5D%29%29%20AND%20%28identity%3A0.9%29"
    url2 = "https://rest.uniprot.org/uniref/stream?compressed=true&fields=id%2Cname%2Ctypes%2Ccount&format=tsv&query=%28%28taxonomy_id%3A2759%29%20AND%20%28identity%3A0.9%29%20AND%20%28count%3A%5B3%20TO%20%2A%5D%29%29""
    request_api(url1, "./data/uniref/uniref_ids_1.tsv")
    request_api(url2, "./data/uniref/uniref_ids_2.tsv")
    merge_outputs("./data/uniref/uniref_ids_1.tsv", "./data/uniref/uniref_ids_2.tsv", "./data/uniref/uniref_ids.tsv")