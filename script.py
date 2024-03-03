import re

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


REFERENCE_LABEL = "NC_045512.2"
genomes = dict()
f = open("aligned-sequences.fasta", 'r')
index = ""
for line in f:
    if line[0] == ">":
        index = re.search(r'>(\S+)', line).group(1)
        genomes[index] = ""
    else:
        genomes[index] += line[:-1]

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

for key in variants_dict:
    print(key)
    print(variants_dict[key])

(a, b) = get_significant_genomes(variants_dict)
array = get_common_mutations(a, variants_dict, True)
print(array)



