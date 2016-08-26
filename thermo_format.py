#!/usr/bin/env python3

# Import modules
import argparse

# Define functions
def format_thermo_lines(thermo_lines):

    # Create dictionary by metabolite
    thermo_lines_dict = {}
    for line in thermo_lines:
        line = line.strip()
        if line.endswith("nan"):
            # Skip lines without a formation delta G
            continue
        if line.startswith("Compound ID"):
            # Skip header line
            continue
        metabolite = line.split(",")[0]
        try:
            thermo_lines_dict[metabolite].append(line)
        except KeyError:
            thermo_lines_dict[metabolite] = [line]

    # Construct return list entry for each metabolite
    net_thermo_lines = []
    for metabolite in sorted(thermo_lines_dict):
        met_section = [metabolite + ";;;;;;"]
        for met_line in thermo_lines_dict[metabolite]:
            met_line = met_line.split(",")
            dfG = met_line[3]
            chrg = met_line[2]
            nH = met_line[1]
            met_section.append(";".join(["", dfG, "NaN", chrg, nH, "", ""]))
        net_thermo_lines.append("\n".join(met_section) + "\n")

    return net_thermo_lines

def test_format_thermo_lines():
    input_thermo_lines = [
        "Compound ID,nH,charge,dG0_f\n",
        "C00008,14,-1,-1974.3299999999999\n",
        "C00008,15,0,-1992.5899999999999\n",
        "C00008,16,1,-2002.6800000000001\n",
        "C00009,0,-3,-1020.02\n",
        "C00009,1,-2,-1093.6099999999999\n",
        "C00009,2,-1,-1133.26\n",
        "C00009,3,0,-1143.53\n",
        "C00010,31,-5,-3026.2800000000002\n",
        "C00010,32,-4,-3083.96\n",
        "C00010,33,-3,-3119.6199999999999\n",
        "C00012,0,0,nan\n"
    ]
    output_thermo_lines = [
        "\n".join([
            "C00008;;;;;;",
            ";-1974.3299999999999;NaN;-1;14;;",
            ";-1992.5899999999999;NaN;0;15;;",
            ";-2002.6800000000001;NaN;1;16;;"
        ]) + "\n",
        "\n".join([
            "C00009;;;;;;",
            ";-1020.02;NaN;-3;0;;",
            ";-1093.6099999999999;NaN;-2;1;;",
            ";-1133.26;NaN;-1;2;;",
            ";-1143.53;NaN;0;3;;"
        ]) + "\n",
        "\n".join([
            "C00010;;;;;;",
            ";-3026.2800000000002;NaN;-5;31;;",
            ";-3083.96;NaN;-4;32;;",
            ";-3119.6199999999999;NaN;-3;33;;"
        ]) + "\n"
    ]
    assert format_thermo_lines(input_thermo_lines) == output_thermo_lines


# Main code block

def main(cc_thermo_filename, net_thermo_filename):
    with open(net_thermo_filename, 'w') as outfile:
        outfile.write(
            "".join(format_thermo_lines(open(cc_thermo_filename).readlines()))
        )

if __name__ == "__main__":

    # Read arguments from the commandline
    parser = argparse.ArgumentParser()

    # Required input: Component-contribution csv file
    parser.add_argument(
        'infile',
        help='Read component-contribution csv file with dfG values.'
    )

    # Output: NET-compatible metabolite thermodynamics file
    parser.add_argument(
        'outfile',
        help='Write NET-formatted metabolite thermodynamics file.'
    )

    args = parser.parse_args()

    # Run main function
    main(args.infile, args.outfile)
