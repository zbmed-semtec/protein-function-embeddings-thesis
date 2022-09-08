from rdflib import Graph, RDF
from rdflib.namespace import RDFS, OWL
from rdflib import Namespace


def create_dictionary(input_file: str) -> dict:
    """Converts input XRDF file into a dictionary with IDs as keys and
    label, aspect, and synonyms as values.
    Parameters
    ----------
    input_file : str
        File path for the input Gene Ontology XRDF file.
    Returns
    ----------
    dic : dict
        Dictionary with IDs as keys and labels, synonyms as values.
    """
    graph = Graph()
    graph.parse(input_file, format="xml")

    go_metadata_dict = {}
    main_namespace = Namespace("http://www.w3.org/2002/07/owl#")
    OBO_OWL = Namespace("http://www.geneontology.org/formats/oboInOwl#")

    for Class in graph.subjects(RDF.type, OWL.Class):
        if Class.startswith("http://purl.obolibrary.org/obo/GO_"):
            deprecated_value = False
            for deprecated in graph.objects(Class, main_namespace.deprecated):
                deprecated_value = deprecated
            if not deprecated_value:
                go_metadata_dict[str(Class)] = {}
                for label in graph.objects(Class, RDFS.label):
                    go_metadata_dict[str(Class)]["label"] = str(label)
                go_metadata_dict[str(Class)]["synonyms"] = []
                go_metadata_dict[str(Class)]["aspect"] = ""
                for synonyms in graph.objects(Class, OBO_OWL.hasExactSynonym):
                    go_metadata_dict[str(Class)]["synonyms"].append(str(synonyms))
                for aspect in graph.objects(Class, OBO_OWL.hasOBONamespace):
                    go_metadata_dict[str(Class)]["aspect"] += str(aspect)
    return go_metadata_dict


def replace_chars(text) -> str:
    """Replaces special characters that invalidate the mwt format with the correct syntax.
    Parameters
    ----------
    text : str
        Term and the Synonyms.
    """
    text = text.replace("&", "&amp;")
    text = text.replace('"', "&quot;")
    text = text.replace("'", "&apos;")
    text = text.replace("<", "&lt;")
    text = text.replace(">", "&gt;")
    return text


def create_mwt_file(dictionary: dict, namespace: str, output_file: str) -> None:
    """Creates a MWT file from a dictionary with IDs as keys and labels, synonyms as values.
    saves the file to the output_file destination.
    Parameters
    ----------
    dictionary: dict
        Gene ontology metadata dictionary.
    namespace: str
        Namespace for the output mwt file.
    output_file: str
        File path for the ouput mwt file.
    """
    with open(output_file, 'w') as output:
        output.write("<?xml version='1.0' encoding='UTF-8'?>\n")
        output.write('<mwt xmlns:z="{}">\n'.format(namespace))
        output.write("<template><z:GO id='%1' aspect='%2'>%0</z:GO></template>\n\n")

        for id, metadata in dictionary.items():
            if 'aspect' in metadata:
                # aspect = ", ".join(map(str, metadata["aspect"]))
                metadata['label'] = replace_chars(metadata["label"])
                output.write('<t p1="{}" p2="{}">{}</t>\n'.
                             format(id, metadata["aspect"], metadata['label']))
                if 'synonyms' in metadata and len(metadata["synonyms"]) != 0:
                    n = 0
                    while n < len(metadata['synonyms']):
                        synonym = metadata['synonyms'][n]
                        synonym = replace_chars(synonym)
                        output.write('<t p1="{}" p2="{}">{}</t>\n'.format(id, metadata["aspect"], synonym))
                        n = n + 1
        output.write("\n</mwt>")
    return


if __name__ == "__main__":
    go_dict = create_dictionary("./data/input/owlapi.xrdf")
    create_mwt_file(go_dict, "https://github.com/zbmed-semtec/protein-function-embeddings-thesis#", "./data/output/go.mwt")
