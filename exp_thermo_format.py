#!/usr/bin/env python3

# Import modules
import argparse
import re

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
    net_thermo_lines = {}
    for metabolite in thermo_lines_dict:
        met_section = [metabolite + ";;;;;;"]
        for met_line in thermo_lines_dict[metabolite]:
            met_line = met_line.split(",")
            dfG = met_line[3]
            chrg = met_line[2]
            nH = met_line[1]
            met_section.append(";".join(["", dfG, "NaN", chrg, nH, "", ""]))
        net_thermo_lines[metabolite] = "\n".join(met_section) + "\n"

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
    output_thermo_lines = {
        "C00008" : "\n".join([
            "C00008;;;;;;",
            ";-1974.3299999999999;NaN;-1;14;;",
            ";-1992.5899999999999;NaN;0;15;;",
            ";-2002.6800000000001;NaN;1;16;;"
        ]) + "\n",
        "C00009" : "\n".join([
            "C00009;;;;;;",
            ";-1020.02;NaN;-3;0;;",
            ";-1093.6099999999999;NaN;-2;1;;",
            ";-1133.26;NaN;-1;2;;",
            ";-1143.53;NaN;0;3;;"
        ]) + "\n",
        "C00010" : "\n".join([
            "C00010;;;;;;",
            ";-3026.2800000000002;NaN;-5;31;;",
            ";-3083.96;NaN;-4;32;;",
            ";-3119.6199999999999;NaN;-3;33;;"
        ]) + "\n"
    }
    assert format_thermo_lines(input_thermo_lines) == output_thermo_lines


def net_model_to_metabolite_dict(net_model_text):
    metabolite_dict = {}
    for line in net_model_text.split("\n"):
        if line.startswith("reaction"):
            equation = line.split(";")[2]
            all_in_one = re.search('^\[.+?\]', equation)
            if all_in_one:
                compartment = all_in_one.group()
            for metabolite in [
                x.strip() for x in re.findall('C[0-9]{5}(?:(?:\[?.+?\]?)|(?= )|$)', equation)
                ]:
                if not all_in_one:
                    compartment = re.search('\[.+\]', metabolite).group()
                    metabolite = metabolite.split("[")[0]
                try:
                    metabolite_dict[metabolite].add(compartment)
                except KeyError:
                    metabolite_dict[metabolite] = {compartment}
    return metabolite_dict

def test_net_model_to_metabolite_dict():
    net_model_text = "\n".join([
        ";ID;pH;IS;Potential mV;Volume;",
        "compartment;c;7.4;0.1;0;0.709;cytosol",
        "compartment;a;7.4;0.1;0;0.007;carboxysome",
        "compartment;e;7.8;0.1;0;0;extracellular",
        "compartment;l;5.0;0.1;0;0.09;thylakoid lumen",
        "compartment;t;7.4;0.1;0;0.06;thylakoid membrane",
        "compartment;y;7.4;0.1;0;0.014;plasma membrane",
        "compartment;p;7.8;0.1;0;0.12;periplasmic space",
        "",
        ";Model;;;;;",
        "reaction;R47;[c](2) C00002 + C00064 + C00288 + C00001 + C00080 = (2) C00008 + C00009 + C00025 + C00169;;;;",
        "reaction;R644;[c]C03561 + C00006 = C00332 + C00005 + C00080;;;;",
        "reaction;R705;C00076[e] + C00002[c] + C00001[c] = C00076[c] + C00008[c] + C00009[c];;;;",
        "reaction;R139;[c]C11355 = C00568 + C00022;;;;",
        "reaction;R48;[c]C00002 + C00014 + C00011 = C00008 + C00169;;;;",
        "reaction;R661;[c]C00353 + (3) C00080 + (3) C00005 = C05427 + (3) C00006;;;;",
        "reaction;R320;[c]C05752 + C01209 = C05753 + C00011 + C00229;;;;",
        "reaction;R322;[c]C05755 + C01209 = C05756 + C00011 + C00229;;;;",
        "reaction;R446;C00005[c] + (5) C00080[c] + C10385[t] = C00006[c] + C16695[t] + (4) C00080[l];;;;",
        "reaction;R544;[c]C00130 + C00001 = C04734;;;;",
        "reaction;R660;[c]C00448 + (8) C00129 = C04574 + (8) C00013;;;;",
        "reaction;R401;[c]C00064 + C00085 = C00025 + C00352;;;;",
        "reaction;R99;[c]C02094 + C00004 + C00080 + C00007 = C08591 + C00003 + C00001;;;;",
        "reaction;R515;[c]C00158 + C00010 = C00024 + C00001 + C00036;;;;",
        "reaction;R743;[c]C06319 + C00006 = C06320 + C00005 + C00080;;;;",
        "",
        "",
        "Thermo names;;",
        "",
        ";Metabolite (don't change);Name in model",
        "metabolite;C00001;C00001",
        "metabolite;C00002;C00002",
        "metabolite;C00003;C00003",
        "metabolite;C00004;C00004",
        "metabolite;C00005;C00005",
        "metabolite;C00006;C00006",
        "metabolite;C00007;C00007",
        "metabolite;C00008;C00008",
        "metabolite;C00009;C00009",
        "metabolite;C00010;C00010",
        "metabolite;C00011;C00011",
        "metabolite;C00013;C00013",
        "metabolite;C00014;C00014",
        "metabolite;C00015;C00015",
        "metabolite;C00016;C00016"
    ]) + "\n"
    metabolite_dict = {
        'C00001' : {'[c]'}, 'C00002' : {'[c]'}, 'C00003' : {'[c]'},
        'C00004' : {'[c]'},
        'C00005' : {'[c]'}, 'C00006' : {'[c]'}, 'C00007' : {'[c]'},
        'C00008' : {'[c]'}, 'C00009' : {'[c]'}, 'C00010' : {'[c]'},
        'C00011' : {'[c]'}, 'C00013' : {'[c]'}, 'C00014' : {'[c]'},
        'C00022' : {'[c]'}, 'C00024' : {'[c]'}, 'C00025' : {'[c]'},
        'C00036' : {'[c]'}, 'C00064' : {'[c]'},
        'C00076' : {'[c]', '[e]'},
        'C00080' : {'[c]', '[l]'},
        'C00085' : {'[c]'}, 'C00129' : {'[c]'}, 'C00130' : {'[c]'},
        'C00158' : {'[c]'}, 'C00169' : {'[c]'}, 'C00229' : {'[c]'},
        'C00288' : {'[c]'}, 'C00332' : {'[c]'}, 'C00352' : {'[c]'},
        'C00353' : {'[c]'}, 'C00448' : {'[c]'}, 'C00568' : {'[c]'},
        'C01209' : {'[c]'}, 'C02094' : {'[c]'}, 'C03561' : {'[c]'},
        'C04574' : {'[c]'}, 'C04734' : {'[c]'}, 'C05427' : {'[c]'},
        'C05752' : {'[c]'}, 'C05753' : {'[c]'}, 'C05755' : {'[c]'},
        'C05756' : {'[c]'}, 'C06319' : {'[c]'}, 'C06320' : {'[c]'},
        'C08591' : {'[c]'}, 'C10385' : {'[t]'}, 'C11355' : {'[c]'},
        'C16695' : {'[t]'}
    }
    assert metabolite_dict == net_model_to_metabolite_dict(net_model_text)


def format_concentrations(conc_text, net_model_dict):
    formatted_lines = [";Metabolite;Lowest Concentration;Highest Concentration;;;"]
    for line in conc_text.split("\n"):
        if line == "":
            # Skip empty line
            continue
        line = line.split("\t")
        metabolite = line[0]
        lo = str(float(line[1])*1000)
        hi = str(float(line[2])*1000)
        try:
            compartments = sorted(net_model_dict[metabolite])
            for compartment in compartments:
                mc = metabolite + compartment
                formatted_lines.append(";".join([
                    "metabolite", mc, lo, hi, "", "", ""
                ]))
        except KeyError:
            continue
    return "\n".join(formatted_lines) + "\n"

def test_format_concentrations():
    conc_text = "\n".join([
        "C00001\t1\t1",
        "C00002\t3.01207422253299e-05\t0.043431863332143",
        "C00003\t0.000131423040406459\t0.000272201955163885",
        "C00006\t0.000125010602979862\t0.00039712588682272",
        "C00007\t0.000230\t0.000230",
        "C00008\t4.27715338667937e-05\t0.00418670269257234",
        "C00009\t0.017\t0.017",
        "C00011\t0.000010\t0.000010",
        "C00016\t6.15157307144173e-06\t0.000102762615216182",
        "C00019\t1.98544815403941e-05\t0.000972460600392821"
    ]) + "\n"
    net_model_dict = {
        'C00001' : {'[c]', '[a]'}, 'C00002' : {'[c]'},
        'C00003' : {'[c]'}, 'C00004' : {'[c]'},
        'C00005' : {'[c]'}, 'C00006' : {'[c]'},
        'C00007' : {'[c]'}, 'C00008' : {'[c]'},
        'C00009' : {'[c]'}, 'C00010' : {'[c]'},
        'C00011' : {'[c]', '[a]'}, 'C00013' : {'[c]'},
        'C00014' : {'[c]'}, 'C00022' : {'[c]'},
        'C00024' : {'[c]'}, 'C00025' : {'[c]'},
    }
    exp_conc_text = "\n".join([
        ";Metabolite;Lowest Concentration;Highest Concentration;;;",
        "metabolite;C00001[a];1000.0;1000.0;;;",
        "metabolite;C00001[c];1000.0;1000.0;;;",
        "metabolite;C00002[c];0.0301207422253299;43.431863332143;;;",
        "metabolite;C00003[c];0.131423040406459;0.272201955163885;;;",
        "metabolite;C00006[c];0.125010602979862;0.39712588682272;;;",
        "metabolite;C00007[c];0.23;0.23;;;",
        "metabolite;C00008[c];0.0427715338667937;4.18670269257234;;;",
        "metabolite;C00009[c];17.0;17.0;;;",
        "metabolite;C00011[a];0.01;0.01;;;",
        "metabolite;C00011[c];0.01;0.01;;;",
    ]) + "\n"
    assert format_concentrations(conc_text, net_model_dict).split("\n") == \
    exp_conc_text.split("\n")


# Main code block

def main(model, thermo, conc, ratios, fluxes, experimental, out_thermo):

    # Read input files
    ther_dict = format_thermo_lines(open(thermo).readlines())
    meta_dict = net_model_to_metabolite_dict(open(model).read())
    conc_text = format_concentrations(open(conc).read(), meta_dict)
    rats_list = [x.strip() for x in open(ratios).readlines()]
    flux_list = [x.strip().split("\t") for x in open(fluxes).readlines()]

    # Generate output text
    x_out = []
    x_out.append(conc_text.strip())
    x_out.extend(["metabolite;" + x for x in rats_list])
    x_out.append("")
    x_out.append(";ID;direction;")
    x_out.extend([";".join(["flux"] + x + [""]) for x in flux_list])
    x_out = "\n".join(x_out + [""])

    metabolites = filter(lambda x: x in ther_dict, sorted(meta_dict))
    t_out = "".join([ther_dict[m] for m in metabolites])

    # Write output to files
    with open(experimental, 'w') as x:
        x.write(x_out)

    with open(out_thermo, 'w') as t:
        t.write(t_out)


if __name__ == "__main__":

    # Read arguments from the commandline
    parser = argparse.ArgumentParser()

    # Input: Model, thermodynamics data, concentration bounds, ratios, fluxes
    parser.add_argument(
        '-m', '--model',
        help='Read NET model file.'
    )
    parser.add_argument(
        '-t', '--thermo',
        help='Read CC thermo file.'
    )
    parser.add_argument(
        '-c', '--concentrations',
        help='Read concentrations file.'
    )
    parser.add_argument(
        '-r', '--ratios',
        help='Read ratios file.'
    )
    parser.add_argument(
        '-f', '--fluxes',
        help='Read fluxes file.'
    )

    # Output: Experimental data file, model-specific thermodynamics file
    parser.add_argument(
        '-e', '--experimental',
        help='Write experimental data file.'
    )
    parser.add_argument(
        '-o', '--out_thermo',
        help='Write model-specific thermo file.'
    )

    args = parser.parse_args()

    # Run main function
    main(
        args.model, args.thermo, args.concentrations, args.ratios,
        args.fluxes, args.experimental, args.out_thermo
    )
