#!/usr/bin/env python3

# Import modules
import argparse
import re
import os

# Specify path to repository
global repo_dir
repo_dir = os.path.dirname(__file__)

# Define functions
def extract_concentrations(net_output_text, net_kegg_dict, kegg_name_dict, label):
    # Parse metabolite concentration lines
    return_lines = [
        ["Label", "ID", "Name", "Compartment", "KEGGID",
         "LowIn", "HighIn", "LowOut", "HighOut"]
    ]
    store_lines = False
    for line in net_output_text.split("\n"):
        if line.startswith("CONCENTRATIONS"):
            store_lines = True
            continue
        if line.startswith("THERMODYNAMIC DATA"):
            return_lines.append("")
            break
        if store_lines and line and not line.startswith("Metabolite;"):
            # This is a line that we want
            line = line.split(";")
            ID = re.sub("\[.*\]", "", line[0])
            try:
                kegg_id = net_kegg_dict[ID]
            except KeyError:
                try:
                    kegg_name = kegg_name_dict[ID]
                    kegg_id = ID
                except KeyError:
                    continue
            try:
                kegg_name = kegg_name_dict[kegg_id]
            except KeyError:
                kegg_name = "NA"
            try:
                compartment = re.findall("\[.*\]", line[0])[0].strip("[]")
            except IndexError:
                compartment = "NA"
            return_lines.append([
                label, ID, kegg_name, compartment, kegg_id, *line[2:6]
            ])
    return "\n".join(["\t".join(L) for L in return_lines])

def test_extract_concentrations():
    net_output_text = open(os.path.join(repo_dir, "data/example_net.csv")).read()
    kegg_name_dict = dict([
        x.strip().split(";")[0].split("\t") for x in open(os.path.join(repo_dir, \
        "data/keggid_keggname.tab")).readlines()
    ])
    net_kegg_dict = {
        'oaa':'C00036', 'prpp':'C00119', 'sbt6p':'C01096',
        'icit':'C00311', 'xylu-D':'C00310', 'nad':'C00003',
        'lcts':'C00243', 'ara5p':'C01112', '4mop':'C00233',
        'ppi':'C00013', 'g6p':'C00092', 'mnl1p':'C00644',
        'no2':'C00088', '6pgl':'C01236', 'ade':'C00147'
    }
    exp_extract = "\n".join([
        "\t".join(["Label", "ID",     "Name",                                    "Compartment", "KEGGID", "LowIn", "HighIn", "LowOut", "HighOut"]),
        "\t".join(["Test",  '4mop',   '4-Methyl-2-oxopentanoate',                'c',           'C00233', "0.0001", "10", "0.0001", "10"]),
        "\t".join(["Test",  '6pgl',   'D-Glucono-1,5-lactone 6-phosphate',       'c',           'C01236', "0.567", "0.693", "0.567", "0.693"]),
        "\t".join(["Test",  'ade',    'Adenine',                                 'c',           'C00147', "0.0001", "10", "0.0001", "10"]),
        "\t".join(["Test",  'ara5p',  'D-Arabinose 5-phosphate',                 'c',           'C01112', "0.0001", "10", "0.0001", "0.253032"]),
        "\t".join(["Test",  'g6p',    'D-Glucose 6-phosphate',                   'c',           'C00092', "0.0001", "0.528", "0.337486", "0.5279"]),
        "\t".join(["Test",  'icit',   'Isocitrate',                              'c',           'C00311', "0.0001", "10", "0.0001", "0.684001"]),
        "\t".join(["Test",  'lcts',   'Lactose',                                 'c',           'C00243', "0.0001", "10", "0.0001", "10"]),
        "\t".join(["Test",  'mnl1p',  'D-Mannitol 1-phosphate',                  'c',           'C00644', "0.0001", "10", "0.0001", "10"]),
        "\t".join(["Test",  'nad',    'NAD+',                                    'c',           'C00003', "0.0001", "10", "0.0331983", "10"]),
        "\t".join(["Test",  'no2',    'Nitrite',                                 'c',           'C00088', "0.0001", "10", "0.0001", "10"]),
        "\t".join(["Test",  'oaa',    'Oxaloacetate',                            'c',           'C00036', "0.0001", "10", "0.00051201", "10"]),
        "\t".join(["Test",  'ppi',    'Diphosphate',                             'c',           'C00013', "1", "100", "1", "100"]),
        "\t".join(["Test",  'prpp',   '5-Phospho-alpha-D-ribose 1-diphosphate',  'c',           'C00119', "0.0001", "10", "0.0001", "10"]),
        "\t".join(["Test",  'sbt6p',  'Sorbitol 6-phosphate',                    'c',           'C01096', "0.0001", "10", "0.0001", "10"]),
        "\t".join(["Test",  'xylu-D', 'D-Xylulose',                              'c',           'C00310', "0.0001", "10", "0.0001", "10"]),
        ""
    ])
    X = extract_concentrations(net_output_text, net_kegg_dict, kegg_name_dict, "Test")
    assert X.split("\n") == exp_extract.split("\n")

# Main code block

def main(infile, names, label, outfile_name):
    kegg_name_dict = dict([
        x.strip().split(";")[0].split("\t") for x in \
        open(os.path.join(repo_dir, "data/keggid_keggname.tab")).readlines()
    ])
    if names:
        net_kegg_dict = dict([
            x.strip().split("\t") for x in open(names).readlines()
        ])
    else:
        net_kegg_dict = {}
    net_output_text = open(infile).read()
    with open(outfile_name, 'w') as outfile:
        output = extract_concentrations(
            net_output_text, net_kegg_dict, kegg_name_dict, label
        )
        outfile.write(output)

if __name__ == "__main__":

    # Read arguments from the commandline
    parser = argparse.ArgumentParser()

    # Input: NET output file, metabolite name translation table
    parser.add_argument(
        '-i', '--infile',
        help='Read NET output file.'
    )
    parser.add_argument(
        '-n', '--names',
        help='Read metabolite names-to-KEGG ID translation.'
    )
    parser.add_argument(
        '-l', '--label',
        help='Label for dataset.'
    )

    # Output: Tab-delimited file with concentration ranges
    parser.add_argument(
        '-o', '--outfile',
        help='Write concentration ranges to outfile.'
    )

    args = parser.parse_args()

    # Run main function
    main(args.infile, args.names, args.label, args.outfile)
