import xml
import pandas as pd
from xml.etree.ElementTree import Element, SubElement
from xml.etree.ElementTree import ElementTree


def write_xml(input_tsv_file: str, output_xml_file: str) -> None:
    """
    Creates and writes an XML file with each Uniprot-SwissProt protein as an entry along with its corresponding metadata.
    Parameters
    ----------
    input_tsv_file: str
        File path for TSV file with protein entries and metadata.
    output_xml_file: str
        File path for output XML file.
    """
    root = Element('collection')
    df = pd.read_csv(input_tsv_file, delimiter='\t', index_col=0, converters={'taxon_lineage': pd.eval})

    for index, row in df.iterrows():
        entries = SubElement(root, 'entry')
        accession = SubElement(entries, 'accession')
        accession.text = str(row['accession'])
        function = SubElement(entries, 'function')
        function.text = str(row['function'])
        evidence = SubElement(entries, 'evidence')
        evidence_tags = str(row['evidence_tags']).split(', ')
        for tag in evidence_tags:
            evidence_tag = SubElement(evidence, 'evidence_tag')
            evidence_tag.text = tag
        evidence_text = SubElement(evidence, 'evidence_text')
        evidence_text.text = str(row['evidence_text'])
        taxonomy = SubElement(entries, 'taxonomy')
        taxon = SubElement(taxonomy, 'taxon')
        taxon.text = str(row['taxon'])
        lineage_list = row['taxonomy_lineage'].strip('[').strip(']').split(', ')
        for lineage in lineage_list:
            taxonomy_lineage = SubElement(taxonomy, 'lineage')
            taxonomy_lineage.text = lineage.strip("'")
        taxon_id = SubElement(taxonomy, 'taxon_id')
        taxon_id.text = str(row['taxon_identifier'].split("NCBI Taxonomy ")[1])

    et = ElementTree(root)
    xml.etree.ElementTree.indent(et, space=" ", level=0)
    with open(output_xml_file, 'wb') as f:
        et.write(f, encoding="utf-8", xml_declaration=True)


if __name__ == "__main__":
    write_xml("./data/output/functions/rev-20220525-UniProtKB.tsv", "./data/output/functions/swissprot.xml")