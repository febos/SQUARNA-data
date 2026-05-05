

res = []

for dataset, tool, M, B in [("S01",  "SQUARNA", m/10, -b/10) for m in range(13,24) for b in range(1,12)]:

    try:
        with open("{}_{}_{}_{}.tsv".format(dataset,tool,M,B)) as file:

            fs = []

            for line in file.readlines()[1:]:

                fs.append(float(line.split('\t')[9]))

        res.append((M,B,round(sum(fs)/len(fs),4)))
    except:
        continue


for r in sorted(res,key = lambda x: -x[2], reverse = True):
    print(*r)
