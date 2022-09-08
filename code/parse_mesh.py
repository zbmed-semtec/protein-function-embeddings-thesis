from rdflib import Graph, RDF
from rdflib.namespace import SKOS, OWL
from rdflib import Namespace


def create_dictionary(input_file: str) -> dict:
    """Converts input TTL file into a dictionary with IDs as keys and
    label, synonyms, semantic types and cui as values.
    Parameters
    ----------
    input_file : str
        File path for the input MeSH TTL file.
    Returns
    ----------
    dic : dict
        Dictionary with IDs as keys and labels, synonyms, semantic types and cui as values.
    """
    graph = Graph()
    graph.parse(input_file, format="turtle")
    mesh_metadata_dict = {}
    umls = Namespace("http://bioportal.bioontology.org/ontologies/umls/")
    for owlClass in graph.subjects(RDF.type, OWL.Class):
        if owlClass.startswith("http://purl.bioontology.org/ontology/MESH/"):
            for notation in graph.objects(owlClass, SKOS.notation):
                mesh_metadata_dict[str(owlClass)] = {}
            for label in graph.objects(owlClass, SKOS.prefLabel):
                mesh_metadata_dict[str(owlClass)]["label"] = str(label)
            mesh_metadata_dict[str(owlClass)]["synonyms"] = []
            mesh_metadata_dict[str(owlClass)]["semantic_types"] = []
            mesh_metadata_dict[str(owlClass)]["cui"] = []
            for synonyms in graph.objects(owlClass, SKOS.altLabel):
                mesh_metadata_dict[str(owlClass)]["synonyms"].append(str(synonyms))
            for semantic_types in graph.objects(owlClass, umls.hasSTY):
                mesh_metadata_dict[str(owlClass)]["semantic_types"].append(str(semantic_types))
            for cui in graph.objects(owlClass, umls.cui):
                mesh_metadata_dict[str(owlClass)]["cui"].append(str(cui))
    return mesh_metadata_dict


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
        MeSH metadata dictionary.
    namespace: str
        Namespace for the output mwt file.
    output_file: str
        File path for the ouput mwt file.
    """
    with open(output_file, 'w') as output:
        output.write("<?xml version='1.0' encoding='UTF-8'?>\n")
        output.write('<mwt xmlns:z="{}">\n'.format(namespace))
        n_parameters = len(dictionary[max(dictionary, key=lambda v: len(dictionary[v]))])
        if n_parameters > 2:
            output.write("<template><z:MESH id='%1' cui='%2' sty='%3'>%0</z:MESH></template>\n\n")
        else:
            output.write("<template><z:MESH id='%1'>%0</z:MESH></template>\n\n")
            
        for id, metadata in dictionary.items():
            if 'cui' and 'semantic_types' in metadata:
                cui = (", ".join(map(str, metadata["cui"]))).replace("'", "")
                semantic_types = ", ".join(map(str, metadata["semantic_types"]))
                metadata['label'] = replace_chars(metadata["label"])
                output.write('<t p1="{}" p2="{}" p3="{}">{}</t>\n'.format(id, cui, semantic_types, metadata['label']))
                if 'synonyms' in metadata:
                    n = 0
                    while n < len(metadata['synonyms']):
                        synonym = metadata['synonyms'][n]
                        synonym = replace_chars(synonym)
                        output.write('<t p1="{}" p2="{}" p3="{}">{}</t>\n'.format(id, cui, semantic_types, synonym))
                        n = n + 1
            elif "synonyms" in metadata:
                metadata['label'] = replace_chars(metadata['label'])
                output.write('<t p1="{}">{}</t>\n'.format(id, metadata['label']))
                n = 0
                while n < len(metadata['synonyms']):
                    synonym = metadata['synonyms'][n]
                    synonym = replace_chars(synonym)
                    output.write('<t p1="{}">{}</t>\n'.format(id, synonym))
                    n = n + 1
            else:
                metadata['label'] = replace_chars(metadata['label'])
                output.write('<t p1="{}">{}</t>\n'.format(id, metadata['label']))
        output.write("\n</mwt>")
    return
                

if __name__ == "__main__":
    mesh_dict = create_dictionary("./data/input/MESH.ttl")
    create_mwt_file(mesh_dict, "https://github.com/zbmed-semtec/protein-function-embeddings-thesis#", "./data/output/mesh.mwt")
