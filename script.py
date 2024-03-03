import re
import argparse
import sys

def gap_length(sequence, start):
    tmp = None
    result = 0
    if start:
        tmp = re.search(r'^-+', sequence)
    else:
        tmp = re.search(r'-+$', sequence)
    if tmp != None:
        result = len(tmp.group())
    return result

#funzione che ritorna il genoma con più variazioni e quello con meno variazioni rispetto al reference
def get_significant_genomes(variants_dict):
    min = 0
    min_name = ""
    max = 0
    max_name = ""
    for key in variants_dict:
        if len(variants_dict[key]) > max or max_name == "":
            max_name = key
            max = len(variants_dict[key])
        if len(variants_dict[key]) < min or min_name == "":
            min_name = key
            min = len(variants_dict[key])
    return (min_name, max_name)

#funzione che ritorna le posizioni del reference rispe]o a cui tutti gli altri genomi variano 
#same = True --> se si vogliono variazioni comuni uguali, same = False --> se si vogliono variazioni comuni
def get_common_mutations(genome_label, variants_dict, same):
    valid = True
    result = list()
    for index in variants_dict[genome_label]:
        valid = True
        for genome in variants_dict:
           # if (not same and (index not in variants_dict[genome])) or\
            #    (same and (index not in variants_dict[genome] or\
             #       (variants_dict[genome_label][index] != variants_dict[genome][index]))):
            if (index not in variants_dict[genome]) or\
                (same and (variants_dict[genome_label][index] != variants_dict[genome][index])):
                valid = False
                break
        if valid:
            result.append(index)
    return result

def create_genome_dict():
    genomes = dict()
    file = open("aligned-sequences.fasta", 'r')
    index = ""
    for line in file:
        if line[0] == ">":
            index = re.search(r'>(\S+)', line).group(1)
            genomes[index] = ""
        else:
            genomes[index] += line[:-1]
    return genomes

def create_variants_dict(genomes):
    variants_dict = dict()
    opening_gap = gap_length(genomes[REFERENCE_LABEL], True)
    closing_gap = gap_length(genomes[REFERENCE_LABEL], False)
    for key in genomes:
        if key != REFERENCE_LABEL:
            variants_dict[key] = dict()
            opening_gap = max(opening_gap, gap_length(genomes[key], True))
            closing_gap = max(closing_gap, gap_length(genomes[key], False))
            for index in range(opening_gap, len(genomes[key]) - closing_gap):
                if genomes[key][index] != genomes[REFERENCE_LABEL][index]:
                    #verifico che ci sia una cancellazione e la gestisco
                    if genomes[REFERENCE_LABEL][index] == "-":
                        variants_dict[key][index] = ("-", genomes[key][index])
                    #verifico che ci sia un inserimento e lo gestisco
                    elif genomes[key][index] == "-":
                        variants_dict[key][index] = ("+", genomes[REFERENCE_LABEL][index])
                    #verifico che ci sia una modifica e la gestisco
                    elif genomes[key][index] != 'N' and genomes[REFERENCE_LABEL][index] != 'N':
                        variants_dict[key][index] = ("S", (genomes[REFERENCE_LABEL][index], genomes[key][index]))
    return variants_dict

def write_report(output, report_info):
    variants_dict = report_info[0]
    significant_genomes = report_info[1]
    common_mutations = report_info[2]
    common_equal_mutations = report_info[3]

    output.write("Variazioni puntuali rilevate rispetto al reference:\n")
    for key in variants_dict:
        output.write(key + ":\n")
        for position in variants_dict[key]:
            if variants_dict[key][position][0] == "-":
                output.write("• Posizione " + str(position) + ": cancellazione di " + variants_dict[key][position][1] + "\n")
            elif variants_dict[key][position][0] == "+":
                output.write("• Posizione " + str(position) + ": inserimento di " + variants_dict[key][position][1] + "\n")
            else:
                sub_tuple = variants_dict[key][position][1]
                output.write("• Posizione " + str(position) + ": sostituzione " + sub_tuple[0] + " -> " + sub_tuple[1] + "\n")
    
    output.write("\nIl genoma con meno variazioni rispetto al reference è " + significant_genomes[0] + "\n")
    output.write("Il genoma con più variazioni rispetto al reference è " + significant_genomes[1] + "\n")
    
    output.write("\nPosizioni del reference rispetto a cui tutti gli altri genomi variano:\n")
    for position in common_mutations:
        output.write("• " + str(position) + "\n")
    
    output.write("\nPosizioni del reference rispetto a cui tutti gli altri genomi variano allo stesso modo:\n")
    for position in common_equal_mutations:
        output.write("• " + str(position) + "\n")

def write_report_pdf(report_info):
    from reportlab.pdfgen import canvas


REFERENCE_LABEL = "NC_045512.2"

def main(txt, pdf):

    genomes = create_genome_dict()
    variants_dict = create_variants_dict(genomes)
    (least_mutated, most_mutated) = get_significant_genomes(variants_dict)
    common_mutations = get_common_mutations(least_mutated, variants_dict, False)
    common_equal_mutations = get_common_mutations(least_mutated, variants_dict, True)
    report_info = [variants_dict, (least_mutated, most_mutated), common_mutations, common_equal_mutations]

    if txt:
        f = open("report.txt", "w")
        write_report(f, report_info)
        f.close()
    if pdf:
        write_report_pdf(report_info)
    write_report(sys.stdout, report_info)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--txt', action=argparse.BooleanOptionalAction, default=False,
                            help="Per visualizzare il report in un file in formato txt")
    parser.add_argument('--pdf', action=argparse.BooleanOptionalAction, default=False, 
                            help="Per visualizzare il report in un file in formato pdf")
    args = parser.parse_args()
    txt = args.txt
    pdf = args.pdf
    pdf = True
    main(txt, pdf)


