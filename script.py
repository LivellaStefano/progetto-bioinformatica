import re
import argparse
import sys
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors

REFERENCE_LABEL = "NC_045512.2"

# crea il dizionario genomes dove le chiavi sono le label dei genomi, 
# mentre il valore sone le sequenze di basi.
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

# determina la lunghezza del gap iniziale(start = True) o 
# finale(start = False) di sequence
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

# crea il dizionario variants_dict dove le chiavi sono le label 
# dei genomi, mentre i valori sono a loro volta dizionari in cui
# le chiavi sono le posizioni del genoma dove si hanno mutazioni 
# mentre i valori sono il tipo di mutazione, ovvero:
#1) Inserimento, ("+", Base_Inserita)
#2) Eliminazione, ("-", Base_Eliminata)
#3) Cambiamento ("S", (Reference, Genoma))
def create_variants_dict(genomes):
    variants_dict = dict()
    reference_opening_gap = gap_length(genomes[REFERENCE_LABEL], True)
    reference_closing_gap = gap_length(genomes[REFERENCE_LABEL], False)
    for key in genomes:
        if key != REFERENCE_LABEL:
            variants_dict[key] = dict()
            opening_gap = max(reference_opening_gap, gap_length(genomes[key], True)) 
            closing_gap = max(reference_closing_gap, gap_length(genomes[key], False)) 
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

# ritorna il genoma con più variazioni e quello con meno variazioni rispetto al reference
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

# ritorna le posizioni del reference rispetto a cui tutti gli altri genomi 
# variano (tutti i genomi devono variare)
# same = True --> se si vogliono variazioni comuni uguali
# same = False --> se si vogliono variazioni comuni
def get_common_mutations_1(genome_label, variants_dict, same):
    valid = True
    result = list()
    for index in variants_dict[genome_label]:
        valid = True
        for genome in variants_dict:
            # A => same
            # B => index not in variants_dict[genome]
            # C => variants_dict[genome_label][index] != variants_dict[genome][index]
            # (not A and B) or (A and (B or C))
            # (not A and B) or (A and B) or (A and C)
            # (B and (not A or A)) or (A and C)
            # B or (A and C)
            if (index not in variants_dict[genome]) or\
                (same and (variants_dict[genome_label][index] != variants_dict[genome][index])):
                valid = False
                break
        if valid:
            result.append(index)
    return result

# ritorna le posizioni del reference rispetto a cui gli altri genomi 
# variano nello stesso modo (non tutti i genomi devono variare)
def get_common_mutations_2(variants_dict):
    seen = dict()
    for genome in variants_dict:
        for pos in variants_dict[genome]:
            if pos not in seen:
                seen[pos] = (variants_dict[genome][pos], True)
            elif seen[pos] != variants_dict[genome][pos]:
                seen[pos] = (seen[pos][0], False)
    result = list()
    result = [key for key in seen if seen[key][1] == True]
    result.sort()
    return result

def write_report(output, report_info):
    variants_dict = report_info[0]
    significant_genomes = report_info[1]
    common_mutations = report_info[2]
    common_equal_mutations = report_info[3]
    common_equal_mutations_2 = report_info[4]

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
                
    output.write("\nIl genoma con meno variazioni rispetto al reference e' " + significant_genomes[0] + "\n")
    output.write("Il genoma con piu' variazioni rispetto al reference e' " + significant_genomes[1] + "\n")
    
    output.write("\nPosizioni del reference rispetto a cui tutti gli altri genomi variano:\n")
    for position in common_mutations:
        output.write("• " + str(position) + "\n")
    
    output.write("\nPosizioni del reference rispetto a cui tutti gli altri genomi variano allo stesso modo (Versione 1):\n")
    for position in common_equal_mutations:
        output.write("• " + str(position) + "\n")

    output.write("\nPosizioni del reference rispetto a cui tutti gli altri genomi variano allo stesso modo (Versione 2):\n")
    for position in common_equal_mutations_2:
        output.write("• " + str(position) + "\n")

def print_helper(y_pos, data, c):
    c.setFont("Helvetica", 12)
    for pos in data:
        if y_pos < 70:
            c.showPage()
            y_pos = 730
        y_pos -= 15
        c.drawString(70, y_pos, "• " + str(pos))
    return y_pos

def write_report_pdf(report_info):
    c = canvas.Canvas("report.pdf", pagesize=letter)
    c.setFillColor(colors.black)
    c.setFont("Helvetica-Bold", 30)
    c.drawString(50, 700, "Report")
    c.setTitle("Report")

    c.setFillColor(colors.black)
    c.setFont("Helvetica-Bold", 18)
    c.drawString(50, 660, "Variazioni puntuali rilevate rispetto al reference:")

    variants_dict = report_info[0]
    significant_genomes = report_info[1]
    common_mutations = report_info[2]
    common_equal_mutations = report_info[3]
    common_equal_mutations_2 = report_info[4]

    y_pos = 635
    for key in variants_dict:
        c.setFillColor(colors.black)
        c.setFont("Helvetica-Bold", 14)
        c.drawString(50, y_pos, key + ":")
        for position in variants_dict[key]:
            if y_pos < 70:
                c.showPage()
                y_pos = 730
            y_pos -= 15
            c.setFont("Helvetica", 12)
            if variants_dict[key][position][0] == "-":
                c.drawString(70, y_pos, "• Posizione " + str(position) + ": cancellazione di " + variants_dict[key][position][1])
            elif variants_dict[key][position][0] == "+":
                c.drawString(70, y_pos, "• Posizione " + str(position) + ": inserimento di " + variants_dict[key][position][1])
            else:
                sub_tuple = variants_dict[key][position][1]
                c.drawString(70, y_pos, "• Posizione " + str(position) + ": sostituzione " + sub_tuple[0] + " -> " + sub_tuple[1])
        y_pos -= 20

    y_pos -= 10
    c.setFont("Helvetica-Bold", 15.5)
    c.drawString(50, y_pos, "Il genoma con meno variazioni rispetto al reference è " + significant_genomes[0])
    y_pos -= 20
    c.drawString(50, y_pos, "Il genoma con piu' variazioni rispetto al reference è " + significant_genomes[1])

    y_pos -=30
    c.drawString(50, y_pos, "Posizioni del reference rispetto a cui tutti gli altri genomi variano:")
    y_pos = print_helper(y_pos, common_mutations, c)

    y_pos -=30
    c.setFont("Helvetica-Bold", 13)
    c.drawString(50, y_pos, "Posizioni del reference rispetto a cui tutti gli altri genomi variano allo stesso modo (V1):")
    y_pos = print_helper(y_pos, common_equal_mutations, c)

    y_pos -=30
    c.setFont("Helvetica-Bold", 13)
    c.drawString(50, y_pos, "Posizioni del reference rispetto a cui tutti gli altri genomi variano allo stesso modo (V2):")
    y_pos = print_helper(y_pos, common_equal_mutations_2, c)

    c.save()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--txt', action=argparse.BooleanOptionalAction, default=False,
                            help="Per visualizzare il report in un file in formato txt")
    parser.add_argument('--pdf', action=argparse.BooleanOptionalAction, default=False, 
                            help="Per visualizzare il report in un file in formato pdf")
    args = parser.parse_args()
    txt = args.txt
    pdf = args.pdf
    
    genomes = create_genome_dict()
    variants_dict = create_variants_dict(genomes)
    (least_mutated, most_mutated) = get_significant_genomes(variants_dict)
    common_mutations = get_common_mutations_1(least_mutated, variants_dict, False)
    common_equal_mutations = get_common_mutations_1(least_mutated, variants_dict, True)
    common_equal_mutations_2 = get_common_mutations_2(variants_dict)
    report_info = [variants_dict, (least_mutated, most_mutated), common_mutations, common_equal_mutations, common_equal_mutations_2]
    if txt:
        f = open("report.txt", "w")
        write_report(f, report_info)
        f.close()
    if pdf:
        write_report_pdf(report_info)
    write_report(sys.stdout, report_info)

if __name__ == "__main__":
    main()