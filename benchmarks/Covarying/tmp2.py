import os, glob

def ReadStockholm(stkfile):
    """Parses Stockholm format into three lists and two dicts"""

    seqnames = [] # Sequence names
    seqdict  = {} # Sequence dict with name keys and sequence values
    gcnames  = [] # Structure names
    gcdict   = {} # Structure dict with name keys and structure values
    headers  = [] # Headers list

    try:
        file = open(stkfile)
    except:
        # Non-standard encoding found in some
        # of the Rfam families
        file = open(stkfile, encoding="iso8859-15")

    for line in file:
        if line.startswith('#=GC '): # Structure lines

            linesplit = line.strip().split()
            seq = linesplit[-1]
            name = ' '.join(linesplit[1:-1])

            if name not in gcdict:
                gcnames.append(name)
                gcdict[name] = seq
            else:
                gcdict[name] += seq

        elif line.startswith('#'):
            # Header lines
            headers.append(line)

        elif line.startswith('//'):
            pass
        elif not line.strip():
            pass
        else:
            # Sequence lines
            linesplit = line.strip().split()
            seq = linesplit[-1]
            name = ' '.join(linesplit[:-1])

            if name not in seqdict:
                seqnames.append(name)
                seqdict[name] = seq
            else:
                seqdict[name] += seq

    file.close()

    # Put #=GF lines to the end of the headers
    headers1 = [x for x in headers if not x.startswith("#=GF SQ")]
    headers2 = [x for x in headers if x.startswith("#=GF SQ")]
    headers = headers1 + headers2

    return headers, seqnames, seqdict, gcnames, gcdict





fams = [x.split('/')[-1].split('.')[0] for x in sorted(glob.glob("Rfam/*.sto"))]


with open ("Covarying.tsv",'w') as outp:
    outp.write('\t'.join("FAMILY LENGTH DEPTH STRUCTURE TP FP FN PREC REC FSCORE".split())+'\n')
    for fam in fams:
        print(fam)
        headers, seqnames, seqdict, gcnames, gcdict = ReadStockholm("Rfam/{}.sto".format(fam))

        ln = len(seqdict[seqnames[0]])
        dp = len(seqnames)

        for tool in "Rfam CentroidAlifold RNAalifold RscapeNested RscapeTotal SQUARNAs3u".split():
            print(tool)
            sto = "../{}/{}.sto".format(tool,fam)

            os.system("cd tmp; ~/software/rscape_v2.5.6/bin/R-scape --nofigures --rna {} > tmp.tmp".format(sto))

            with open("tmp/tmp.tmp") as inp:
                for line in inp:
                    if line.startswith("# GTp"):
                        nums = line.split()[6:9]
                        TP, PRED, TRUE = int(nums[0]), int(nums[1]), int(nums[2])
                        FP = PRED - TP
                        FN = TRUE - TP
                        break
                else:
                    print("FCKFCKFCK")

            prc = (round(TP / (TP + FP), 3)) if (TP + FP) else 1
            rcl = (round(TP / (TP + FN), 3)) if (TP + FN) else 1
            fsc = (round(2*TP / (2*TP + FP + FN), 3)) if (2*TP + FP + FN) else 1
            outp.write('\t'.join([str(x) for x in [fam, ln, dp, tool, TP, FP, FN, prc, rcl, fsc]])+'\n')


