#!/usr/bin/env python3

# Load modules
import re
import sys
import argparse
from string import ascii_lowercase

# Define functions
def create_compartment_dict(equations):
    # Count compartment occurrences
    compartment_count = {}
    for eq in equations:
        for cm in re.findall("\[.+?\]", eq):
            cm = cm.strip("[]")
            try:
                compartment_count[cm] += 1
            except KeyError:
                compartment_count[cm] = 1

    # Sort by number of occurrences
    cm_count_inv = {}
    for cm in compartment_count:
        try:
            cm_count_inv[compartment_count[cm]].append(cm)
        except KeyError:
            cm_count_inv[compartment_count[cm]] = [cm]

    # Sort alphabetically
    for count in cm_count_inv:
        cm_count_inv[count] = sorted(cm_count_inv[count])

    # Create new tags
    alphabet = list(ascii_lowercase)
    new_tags = {}
    for count in sorted(cm_count_inv.keys(), reverse=True):
        for old_tag in cm_count_inv[count]:
            # Create a new single-letter tag from letters of the old tag or the
            # alphabet, depending on what is already taken
            for new_tag in list(old_tag) + alphabet:
                if new_tag not in new_tags.values():
                    new_tags[old_tag] = new_tag
                    break
                if new_tag == alphabet[-1]:
                    sys.exit("Error: Too many compartments for the alphabet.")

    # Print key-value pairs to inform about new tags
    tags_changed = False
    for tag_pair in new_tags.items():
        if tag_pair[0] != tag_pair[1]:
            tags_changed = True
            break
    if tags_changed:
        print("Compartment tags have been changed:")
        for tag_pair in new_tags.items():
            print("%s --> %s" % tag_pair)

    # Return translation dictionary
    return new_tags

def test_create_compartment_dict():
    # Input 1
    iJO1366_formatted = [
        "udpgal[e]  <=> udpgal[p] ",
        "glu-L[c] + udpLa4o[c]  <=> akg[c] + udpLa4n[c] ",
        "uri[e]  <=> uri[p] ",
        "atp[c] + h2o[c] + taur[p]  -> adp[c] + h[c] + pi[c] + taur[c] ",
        "4 h[c] + sufbcd-4fe4s[c]  -> 4fe4s[c] + sufbcd[c] ",
        "2 dmlz[c]  -> 4r5au[c] + ribflv[c] ",
        "2 o2[c] + q8h2[c]  -> 2 h[c] + 2 o2s[c] + q8[c] ",
        "2omph[c] + 0.5 o2[c]  -> 2ombzl[c] ",
        "3 q8h2[c] + 2 h[p] + no2[p]  -> 3 q8[c] + 2 h2o[p] + nh4[p] ",
    ]

    # Input 2
    knoop_formatted = [
        "PQH_B6_L_[cym] + 1 HB3p_B6_[cym] => 1 C10385_B6_S_[cym] + 1 C00080_[pps] + 1 HB2p_B6_[cym]",
        "C00254_[cyt] => C00166_[cyt] + C00001_[cyt] + C00011_[cyt]",
        "C01269_[cyt] => C00251_[cyt] + C00009_[cyt]",
        "2 C00430_[cyt] => C00931_[cyt] + 2 C00001_[cny]",
        "C03319_[cyt] + C00006_[cyn] <= C00688_[cyt] + C00005_[cyn] + C00080_[cyt]",
        "1.32535 C00093_[cyt] + 1.3327 C05764_[cyt]  => C00416_PG_[cyt] + 2.6507 C00229_[cyt]",
        "C00011_[cyt] => C00011_[ext]"
    ]

    # Expected output 1
    iJO1366_exp_output = {"e":"e", "p":"p", "c":"c"}

    # Expected output 2
    knoop_exp_output = {
        "cyt":"c", "cym":"y", "ext":"e", "pps":"p", "cyn":"n", "cny":"a"
    }

    assert create_compartment_dict(iJO1366_formatted) == iJO1366_exp_output
    assert create_compartment_dict(knoop_formatted) == knoop_exp_output


def reformat_reaction(equation, name_kegg_dict, compartment_dict):

    # Clean up equation
    equation = equation.strip()

    # Split reaction into left- and right-hand sides
    rl, rr = re.split(" +[\=\-\>\<]+ +", equation)

    # Check compartments
    cms = set([x.strip("[]") for x in re.findall("\[.+?\]", equation)])

    # Split on space and iterate through elements of formula
    def stoich(eq_side):
        new_eq_side = []
        for element in [x.split() for x in re.split(" \+ ", eq_side)]:
            if len(element) not in (1,2):
                sys.exit("Error: Equation is badly formatted (%s)" % equation)
            else:
                try:
                    kegg_id = name_kegg_dict[element[-1]]
                except KeyError:
                    # With no corresponding KEGG ID, return empty string
                    return ""
                cm = compartment_dict[re.search("\[.+?\]", element[-1]).group().strip("[]")]
                if len(cms) > 1:
                    kegg_id = kegg_id + "[" + cm + "]"
                if len(element) > 1:
                    new_eq_side.append("(" + element[0] + ")" + " " + kegg_id)
                else:
                    new_eq_side.append(kegg_id)
        return " + ".join(new_eq_side)

    rl = stoich(rl)
    rr = stoich(rr)

    # Empty equation sides mean lack of KEGG ID; reaction should be ignored
    if rl == "" or rr == "":
        return ""

    equation = rl + " = " + rr

    if len(cms) == 1:
        return "[" + compartment_dict[list(cms)[0]] + "]" + equation
    else:
        return equation

def test_reformat_reaction():

    # Input format
    iJO1366_formatted = [
        "udpgal[e]  <=> udpgal[p] ",
        "glu-L[c] + udpLa4o[c]  <=> akg[c] + udpLa4n[c] ",
        "uri[e]  <=> uri[p] ",
        "atp[c] + h2o[c] + taur[p]  -> adp[c] + h[c] + pi[c] + taur[c] ",
        "4 h[c] + sufbcd-4fe4s[c]  -> 4fe4s[c] + sufbcd[c] ",
        "2 dmlz[c]  -> 4r5au[c] + ribflv[c] ",
        "2 o2[c] + q8h2[c]  -> 2 h[c] + 2 o2s[c] + q8[c] ",
        "2omph[c] + 0.5 o2[c]  -> 2ombzl[c] ",
        "3 q8h2[c] + 2 h[p] + no2[p]  -> 3 q8[c] + 2 h2o[p] + nh4[p] ",
    ]

    # Desired format
    net_formatted = [
        "C00052[e] = C00052[p]",
        "[c]C00025 + C16155 = C00026 + C16153",
        "C00299[e] = C00299[p]",
        "C00002[c] + C00001[c] + C00245[p] = C00008[c] + C00080[c] + C00009[c] + C00245[c]",
        "",
        "[c](2) C04332 = C04732 + C00255",
        "[c](2) C00007 + C00390 = (2) C00080 + (2) C00704 + C00399",
        "[c]C05812 + (0.5) C00007 = C99999",
        "(3) C00390[c] + (2) C00080[p] + C00088[p] = (3) C00399[c] + (2) C00001[p] + C01342[p]"
    ]

    # Name to KEGG dictionary
    name_kegg_dict = {
        "udpgal[e]":"C00052", "udpgal[p]":"C00052", "glu-L[c]":"C00025",
        "udpLa4o[c]":"C16155", "akg[c]":"C00026", "udpLa4n[c]":"C16153",
        "uri[e]":"C00299", "uri[p]":"C00299", "atp[c]":"C00002",
        "h2o[c]":"C00001", "taur[p]":"C00245", "adp[c]":"C00008",
        "taur[c]":"C00245", "h[c]":"C00080", "pi[c]":"C00009",
        "dmlz[c]":"C04332", "o2[c]":"C00007", "q8h2[c]":"C00390",
        "4r5au[c]":"C04732", "ribflv[c]":"C00255", "o2s[c]":"C00704",
        "q8[c]":"C00399", "2omph[c]":"C05812", "2ombzl[c]":"C99999",
        "no2[p]":"C00088", "nh4[p]":"C01342", "h[p]":"C00080", "h2o[p]":"C00001"
    }

    # Compartment dictionary
    cm = create_compartment_dict(iJO1366_formatted)

    # Check each case
    for i in range(len(iJO1366_formatted)):
        assert reformat_reaction(iJO1366_formatted[i], name_kegg_dict, cm) == net_formatted[i]

    # Input format
    knoop_formatted = [
        "C00254_[cyt] => C00166_[cyt] + C00001_[cyt] + C00011_[cyt]",
        "C01269_[cyt] => C00251_[cyt] + C00009_[cyt]",
        "2 C00430_[cyt] => C00931_[cyt] + 2 C00001_[cyt]",
        "C03319_[cyt] + C00006_[cyt] <= C00688_[cyt] + C00005_[cyt] + C00080_[cyt]",
        "PQH_B6_L_[cym] + 1 HB3p_B6_[cym] => 1 C10385_B6_S_[cym] + 1 C00080_[pps] + 1 HB2p_B6_[cym]",
        "1.32535 C00093_[cyt] + 1.3327 C05764_[cyt]  => C00416_PG_[cyt] + 2.6507 C00229_[cyt]",
        "C00011_[cyt] => C00011_[ext]"
    ]

    # Desired format
    net_formatted = [
        "[c]C00254 = C00166 + C00001 + C00011",
        "[c]C01269 = C00251 + C00009",
        "[c](2) C00430 = C00931 + (2) C00001",
        "[c]C03319 + C00006 = C00688 + C00005 + C00080",
        "",
        "[c](1.32535) C00093 + (1.3327) C05764 = C00416 + (2.6507) C00229",
        "C00011[c] = C00011[e]"
    ]

    # Name to KEGG dictionary
    name_kegg_dict = {
        "C00254_[cyt]":"C00254", "C00166_[cyt]":"C00166",
        "C00001_[cyt]":"C00001", "C00011_[cyt]":"C00011",
        "C01269_[cyt]":"C01269", "C00251_[cyt]":"C00251",
        "C00009_[cyt]":"C00009", "C00430_[cyt]":"C00430",
        "C00931_[cyt]":"C00931", "C03319_[cyt]":"C03319",
        "C00006_[cyt]":"C00006", "C00688_[cyt]":"C00688",
        "C00005_[cyt]":"C00005", "C00080_[cyt]":"C00080",
        "C10385_B6_S_[cym]":"C10385", "C00080_[pps]":"C00080",
        "C00093_[cyt]":"C00093", "C05764_[cyt]":"C05764",
        "C00416_PG_[cyt]":"C00416", "C00229_[cyt]":"C00229",
        "C00011_[cyt]":"C00011", "C00011_[ext]":"C00011"
    }

    # Compartment dictionary
    cm = create_compartment_dict(knoop_formatted)

    # Check each case
    for i in range(len(knoop_formatted)):
        assert reformat_reaction(knoop_formatted[i], name_kegg_dict, cm) == net_formatted[i]


# Main code block
def main(metabolites, reactions, compartments, biomass, outfile_name):

    # Read metabolite table into dictionary
    name_kegg_dict = dict(
        [L.strip().split("\t") for L in open(metabolites, 'r').readlines()]
    )

    # Remove entries that do not have a valid KEGG ID
    no_kegg_id_names = set()
    for name_kegg_pair in name_kegg_dict.items():
        if not re.match("^C[0-9]{5}$", name_kegg_pair[1]):
            no_kegg_id_names.add(name_kegg_pair[0])
    for n in no_kegg_id_names:
        del(name_kegg_dict[n])

    # Read reactions
    reaction_dict = dict(
        [L.split("\t") for L in open(reactions, 'r').readlines()]
    )

    # Construct compartment dictionary
    cm_dict = create_compartment_dict(reaction_dict.values())

    # Read compartment description file
    compartment_description = [
        L.split("\t") for L in open(compartments, 'r').readlines()
    ]

    # Construct and write output to outfile
    with open(outfile_name, 'w') as f:

        # Write compartment description to outfile
        f.write(";ID;pH;IS;Potential mV;Volume;\n")
        for line in compartment_description:
            line = "compartment;" + cm_dict[line[0]] + ";" + ";".join(line[1:])
            f.write(line)

        # Write "Model" header to outfile
        f.write("\n")
        f.write(";Model;;;;;\n")
        f.write("\n")

        # Write reactions to outfile
        written_reactions = set()

        f.write(";Abbreviation;reactions;;;;\n")
        for reaction_id in sorted(reaction_dict):
            reaction = reformat_reaction(
                reaction_dict[reaction_id], name_kegg_dict, cm_dict
            )
            if reaction in written_reactions:
                # Do not write more than the first of one specific reaction
                continue
            if reaction:
                f.write("reaction;" + reaction_id + ";" + reaction + ";;;;\n")
                written_reactions.add(reaction)

        # Write biomass reaction to outfile (if applicable)
        if biomass:
            biomass_reaction = re.sub("\[.+?\]", "", reformat_reaction(
                open(biomass, 'r').read().strip(), name_kegg_dict, cm_dict
            ))
            f.write("\n")
            f.write(";Biomass Reaction;\n")
            f.write("reaction;Biomass;" + biomass_reaction + "\n")

        # Write "Thermo names header to outfile"
        f.write("\n")
        f.write("\n")
        f.write("Thermo names;;\n")
        f.write("\n")

        # Write metabolite names to outfile
        f.write(";Metabolite (don't change);Name in model\n")
        for kegg_id in sorted(set(name_kegg_dict.values())):
            f.write("metabolite;" + kegg_id + ";" + kegg_id + "\n")
        f.write("\n")

if __name__ == "__main__":

    # Read arguments from the commandline
    parser = argparse.ArgumentParser()

    # Required input: Tables
    parser.add_argument(
        '-m', '--metabolites', required=True,
        help='Read tab-delimited file with metabolite names and IDs.'
    )
    parser.add_argument(
        '-r', '--reactions', required=True,
        help='Read tab-delimited file with reaction IDs and equations.'
    )
    parser.add_argument(
        '-c', '--compartments', required=True,
        help='Read tab-delimited file with compartment properties.'
    )

    # Optional input: Biomass reaction
    parser.add_argument(
        '-b', '--biomass',
        help='File with one line describing the biomass reaction.'
    )

    # Output: NET-compatible model text file
    parser.add_argument(
        '-o', '--outfile', type=str, required=True,
        help='NET-formatted model text file.'
    )

    args = parser.parse_args()

    # Run main function
    main(
        args.metabolites, args.reactions, args.compartments,
        args.biomass, args.outfile
    )
