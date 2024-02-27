import re

genomes = dict()
f = open("covid-sequences.fasta", 'r')
index = ""
for line in f:
    if line[0] == ">":
        index = re.search(r'>(\S+)', line).group(1)
        genomes[index] = ""
    else:
        genomes[index] += line[:-1]
print(genomes["OL700544.1"])





