

res = []

for dataset, tool, M, B in [("Ribonanza2A3",  "RNAfold", m/10, -b/10) for m in range(0,31) for b in range(-20,21)]:

    with open("{}_{}_{}_{}.tsv".format(dataset,tool,M,B)) as file:

        fs = []

        for line in file.readlines()[1:]:

            fs.append(float(line.split('\t')[9]))

    res.append((M,B,round(sum(fs)/len(fs),4)))


for r in sorted(res,key = lambda x: -x[2]):
    print(*r)
