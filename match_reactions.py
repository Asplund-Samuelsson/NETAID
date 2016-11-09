#!/usr/bin/env python3

# Import modules
import argparse
import re
import os

# Specify path to repository
global repo_dir
repo_dir = os.path.dirname(__file__)

# Define functions
def match(reaction1, reaction2):

    # Fail immediately if reactions are empty
    if not reaction1 or not reaction2:
        return (False, 0)

    # Find out if the reaction is restricted to one compartment
    try:
        compartment1 = re.match("^\[.+\]", reaction1).group()
    except AttributeError:
        compartment1 = ""
    try:
        compartment2 = re.match("^\[.+\]", reaction2).group()
    except AttributeError:
        compartment2 = ""

    # Remove leading compartment specification
    reaction1 = re.sub("^\[.+\]", "", reaction1)
    reaction2 = re.sub("^\[.+\]", "", reaction2)

    # Split reactions
    reaction1 = re.split(" +[\=\-\>\<]+ +", re.sub("^\[.+?\]", "", reaction1))
    reaction2 = re.split(" +[\=\-\>\<]+ +", re.sub("^\[.+?\]", "", reaction2))

    # Format the reactions in a consistent manner
    def format_reaction(reaction, compartment):
        for i in [0,1]:
            reaction[i] = reaction[i].split(" + ")
            for j in range(len(reaction[i])):
                if "(" in reaction[i][j]:
                    reaction[i][j] = re.sub("[()]", "", reaction[i][j])
                    reaction[i][j] = reaction[i][j].split(" ")
                    reaction[i][j][0] = float(reaction[i][j][0])
                    reaction[i][j][1] = reaction[i][j][1] + compartment
                else:
                    reaction[i][j] = [1.0, reaction[i][j] + compartment]
            # Sort the elements by compound
            reaction[i] = sorted(reaction[i], key = lambda x: x[1])
        return reaction

    r1_unsorted = format_reaction(reaction1, compartment1)
    r2_unsorted = format_reaction(reaction2, compartment2)

    # Sort the sides
    reaction1 = sorted(r1_unsorted, key = lambda x: "".join([y[1] for y in x]))
    reaction2 = sorted(r2_unsorted, key = lambda x: "".join([y[1] for y in x]))

    # Get metabolites for comparison
    reaction1_met = [[x[1] for x in y] for y in reaction1]
    reaction2_met = [[x[1] for x in y] for y in reaction2]

    # Get stoichiometry for comparison
    reaction1_sto = [[x[0] for x in y] for y in reaction1]
    reaction2_sto = [[x[0] for x in y] for y in reaction2]

    # Compare reactions and return verdict
    if reaction1_met != reaction2_met:
        return (False, 0)

    sto_set = set()

    for i in [0,1]:
        for j in range(len(reaction1_sto[i])):
            sto_set.add(reaction1_sto[i][j]/reaction2_sto[i][j])
            if len(sto_set) > 1:
                return (False, 0)

    # Reactions are matched; Determine directionality
    r1_met = [[x[1] for x in y] for y in r1_unsorted]
    r2_met = [[x[1] for x in y] for y in r2_unsorted]

    # Reactions matched; same direction
    if r1_met == r2_met:
        return (True, 1)
    # Reactions matched; different direction
    else:
        return (True, -1)


def test_match():
    reactions_1 = [
        "C00084[e] = C00084[c]",
        "[c]C00049 + C00044 + C00130 = C03794 + C00035 + C00009",
        "[c]C01268 + C00005 = C04454 + C00006",
        "[c](2) C00027 = (2) C00001 + C00007",
        "[c]C03263 + C00007 = (2) C00011 + (2) C00001 + C01079",
        "(0.5) C00007[c] + C00390[c] = C00001[c] + C00399[c]",
        "[c]C00026 + C00064 + C00005 = (2) C00025 + C00006",
        "[c]C00311 = C00048 + C00042",
        "[c]C00143 + C00004 = C00440 + C00003",
        "[c]C00092 = C00085",
        "[c]C00399 + C00042 = C00390 + C00122",
        "[c]C00149 + C00003 = C00011 + C00004 + C00022"
    ]
    reactions_2 = [
        "C00084[c] = C00084[e]",
        "[c]C00044 + C00130 + C00049 = C00035 + C00009 + C03794",
        "[c]C04454 + C00006 = C01268 + C00005",
        "[c](2) C00027 = C00007 + (2) C00001",
        "[c]C03263 + C00007 = C01079 + (2) C00011 + (2) C00001",
        "C00007[c] + (2) C00390[c] = (2) C00001[c] + (2) C00399[c]",
        "[c](2) C00025 + C00003 = C00064 + C00026 + C00004",
        "[c]C00311 = C00042 + C00048",
        "[c]C00440 + C00006 = C00143 + C00005",
        "[c]C00092 = C00085",
        "",
        "[c]C00149 + C00006 = C00022 + C00011 + C00005"
    ]
    expected_matches = [
        (True, -1), (True, 1), (True, -1), (True, 1),
        (True, 1), (True, 1), (False, 0), (True, 1),
        (False, 0), (True, 1), (False, 0), (False, 0)
    ]

    for pair in zip(reactions_1, reactions_2, expected_matches):
        assert match(pair[0], pair[1]) == pair[2]
        assert match(pair[1], pair[0]) == pair[2]

# Main code block

def main(infile1, infile2, single, outfile_name):
    # Read reactions from infile 1
    reactions_1 = {}
    for line in open(infile1, 'r').readlines():
        if line.startswith("reaction"):
            if single and ";[" not in line:
                continue
            line = line.strip().split(";")
            reactions_1[line[1]] = line[2]

    # Read reactions from infile 2
    reactions_2 = {}
    for line in open(infile2, 'r').readlines():
        if line.startswith("reaction"):
            if single and ";[" not in line:
                continue
            line = line.strip().split(";")
            reactions_2[line[1]] = line[2]

    # Check all reaction combinations
    outfile = open(outfile_name, 'w')
    for rxn_id_1 in reactions_1:
        for rxn_id_2 in reactions_2:
            match_result = match(reactions_1[rxn_id_1], reactions_2[rxn_id_2])
            if match_result[0]:
                outfile.write("\t".join([rxn_id_1, rxn_id_2, str(match_result[1])]) + "\n")
    outfile.close()

if __name__ == "__main__":

    # Read arguments from the commandline
    parser = argparse.ArgumentParser()

    # Input: NET model infiles
    parser.add_argument(
        'infile1',
        help='Read NET model infile 1.'
    )
    parser.add_argument(
        'infile2',
        help='Read NET model infile 2.'
    )

    # Options
    parser.add_argument(
        '-s', '--single', action='store_true',
        help='Single compartment reactions only.'
    )

    # Output: Tab-delimited file specifying matching reaction pairs
    parser.add_argument(
        'outfile',
        help='Write matched reactions to outfile.'
    )

    args = parser.parse_args()

    # Run main function
    main(args.infile1, args.infile2, args.single, args.outfile)
