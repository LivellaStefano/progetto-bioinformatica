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
    print(variants_dict[key])








