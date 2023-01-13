import pandas as pd
import numpy as np
import logging
import xml.etree.ElementTree as ET
import csv


logging.basicConfig(filename='./data/logging/uniref_clusters.log', filemode='a',
                    level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


def check_elem(elem: ET.Element, euk_elems: np.array) -> bool:
    """
    Checks if the element is present in the TSV file of clusters of interest(Eykaryota, 90% identity, >= 2 members,
    only swissprot entries).
    Parameters
    ----------
    elem : ET.Element
        Current element while parsing the XML.
    euk_elems : np.array
        Numpy array of all Uniref cluster ids of interest.
    Returns
    -------
    bool : True/False
        Returns True if the current cluster element belongs to our clusters of interest.
    """
    id = elem.attrib['id']
    if id in euk_elems:
        return True
    else:
        return False


def extract_cluster_metadata(elem: ET.Element) -> dict:
    """
    Extracts the metadata corresponding to the input cluster element. Returns a dictionary consisting of the Uniref
    cluster id, cluster name, members and cluster size.
    Parameters
    ----------
    elem : ET.Element
        Current cluster element.
    Returns
    -------
    cluster_dict : dict
        Dictionary of cluster metadata, consists of id, name, members and member count (cluster size).
    """
    global member_count, cluster_name
    cluster_id = elem.attrib['id']
    cluster_dict = {}
    list_members = []
    for child in elem:
        # print(child)
        if child.tag == '{http://uniprot.org/uniref}name':
            cluster_name = child.text
        if child.tag == '{http://uniprot.org/uniref}representativeMember':
            for sub_child in child:
                if sub_child.tag == '{http://uniprot.org/uniref}dbReference':
                    for sup_sub_child in sub_child:
                        if 'type' in sup_sub_child.attrib and sup_sub_child.attrib['type'] == 'UniProtKB accession':
                            member = sup_sub_child.attrib['value']
                            list_members.append(member)
        if child.tag == '{http://uniprot.org/uniref}member':
            for sub_child in child:
                if sub_child.tag == '{http://uniprot.org/uniref}dbReference':
                    for sup_sub_child in sub_child:
                        if 'type' in sup_sub_child.attrib and sup_sub_child.attrib['type'] == 'UniProtKB accession':
                            member = sup_sub_child.attrib['value']
                            list_members.append(member)
        if 'type' in child.attrib:
            if child.attrib['type'] == 'member count':
                member_count = child.attrib['value']
    cluster_dict[cluster_id] = {'cluster_name': cluster_name, 'n_members': member_count, 'members': list_members}
    return cluster_dict


def write_clusters(cluster_metadata: dict, output_file: str) -> None:
    """
    Write the cluster metadata as a row in the given TSV file.
    Parameters
    ----------
    cluster_metadata : dict
        Dictionary of cluster metadata, consists of id, name, members and member count (cluster size).
    output_file : str
        Filepath for the TSV file consisting of Uniref clusters.
    """
    tsv_columns = ['cluster_id', 'cluster_name', 'members', 'n_members']
    with open(output_file, 'a') as tsv_file:
        writer = csv.DictWriter(tsv_file, fieldnames=tsv_columns, delimiter='\t')
        for key, val in cluster_metadata.items():
            row = {'cluster_id': key}
            row.update(val)
            writer.writerow(row)


def wrapper_function(filtered_uniref_cluster_ids: str, uniref_xml: str, output_file: str) -> None:
    """
    Wrapper function which reads the input files, parses the Uniref 90% XML file, extracts metadata for cluster of
    interest and writes it as a row onto a TSV file.
    Parameters
    ----------
    filtered_uniref_cluster_ids : str
        Filepath of filtered uniref cluster ids as a TSV file.
    uniref_xml : str
        Filepath to Uniref 90% XML file.
    output_file : str
        Filepath to save the TSV of Uniref clusters.
    """
    euk_elems = pd.read_csv(filtered_uniref_cluster_ids, sep='\t')
    euk_elems = np.array(euk_elems['Cluster ID'])

    context = ET.iterparse(uniref_xml, events=("start", "end"))
    context = iter(context)
    event, root = next(context)

    main_counter = 0
    logging.info("START")

    for event, elem in context:
        if event == "end" and elem.tag == "{http://uniprot.org/uniref}entry":
            main_counter += 1
            result = check_elem(elem, euk_elems)
            if result:
                cluster = extract_cluster_metadata(elem)
                write_clusters(cluster, output_file)
            if main_counter % 100000 == 0:
                logging.info(f"Read {main_counter} entries")
            elem.clear()


if __name__ == "__main__":
    wrapper_function("./data/uniref/filtered_uniref_cluster_ids.tsv",
                     "./data/uniref90.xml",
                     "./data/uniref/uniref_clusters.tsv")
