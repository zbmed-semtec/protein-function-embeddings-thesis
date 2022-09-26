import os
import re
from bs4 import BeautifulSoup
import warnings
import logging

logging.basicConfig(filename='./data/logging/annotation.log', filemode='w',
                    level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


def format_input(input_path, formatted_input_path) -> None:
    """Formats the input XML file by unescaping the escaped special characters.
    Replaces the HTML encoding of the special characters to its original form.
    The conversion is as follows:
    &lt; --> <
    &gt: --> >
    &amp; --> &

    Parameters
    ----------
    input_path : str
        Filepath of the input XML file.
    formatted_input_path : str
        Filepath for the formatted XML file.
    """
    logging.info("Formatting input file")
    output = os.path.basename(formatted_input_path)
    os.system(f"sed 's/&lt;/</g;s/&gt;/>/g;s/&amp;/\&/g' {input_path} > {output}")
    logging.info("Formatted input file")


def annotate(formatted_input_path, output_path) -> None:
    """Annotates the formatted input XML files with the dictionary of choice in server file using the whatizit tool.

    Parameters
    ----------
    formatted_input_path : str
        Filepath of the formatted XML file.
    output_path : str
        Filepath for the annotated XML file.
    """
    os.system(f"cat {formatted_input_path} | DistFilter svr=xmlElem > {output_path}")
    logging.info("Annotated file")


def escape_character(contents) -> str:
    """
    Converts the '<' character to its encoded form.
    Replaces the '<' character if it is followed by alphabets or a period other than those
    in the annotation tags to &lt;.

    Parameters
    ----------
    contents : str
        Body of the annotated XML file.
    Returns
    ----------
    contents : str
        Formatted body of the annotated XML file.
    """
    function = re.findall(r'<function>(.*?)</function>', contents, re.DOTALL)[0]
    # contents = contents.replace(function, re.sub(r"<(?!(/*)z:GO)", "&lt;", function))
    contents = contents.replace(function, re.sub(r"<(?!(/*)z:(MESH|go))", "&lt;", function))
    return contents


def format_output(output_path, formatted_output_path) -> None:
    """Formats the annotated XML output into its valid XML format.

    Parameters
    ----------
    output_path : str
        Filepath of the annotated XML file.
    formatted_output_path : str
        Filepath for the formatted annotated XML file.
    """
    with open(output_path, "r") as file_data:
        # Read each line in the file, readlines() returns a list of lines
        contents = file_data.readlines()
        first_line = contents[0]
        # Adding namespace
        contents[1] = contents[1].split(">")[
                              0] + ' xmlns:z="https://github.com/zbmed-semtec/protein-function-embeddings-thesis#">\n'
        # Combine the lines in the list into a string
        contents = "".join(contents)
        formatted_contents = escape_character(contents)
        soup = BeautifulSoup(formatted_contents, "lxml")
    f = open(formatted_output_path, "w")
    f.write(first_line)
    f.write(str(soup.body.next))
    f.close()
    logging.info("Formatted annotated file")


if __name__ == "__main__":
    format_input("./data/output/functions/swissprot.xml", "./data/output/annotations/swissprot_sed.xml")
    annotate("./data/output/annotations/swissprot_sed.xml", "./data/output/annotations/swissprot_go.xml")
    format_output("./data/output/annotations/swissprot_go.xml", "./data/output/annotations/swissprot_go_formatted.xml")
    format_input("./data/output/functions/swissprot_go_formatted.xml", "./data/output/annotations/swissprot_go_sed.xml")
    annotate("./data/output/annotations/swissprot_go_sed.xml", "./data/output/annotations/swissprot_go_mesh.xml")
    format_output("./data/output/annotations/swissprot_go_mesh.xml", "./data/output/annotations/swissprot_go_mesh_formatted.xml")