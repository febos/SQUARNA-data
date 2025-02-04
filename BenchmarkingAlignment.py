import os, time, glob

def PairsToDBN(newpairs, length = 0, returnlevels = False, levellimit = -1):
    """Convert a list of base pairs into a dbn string of the given length"""

    # Initialize the dbn string
    dbn = ['.']*length

    # Define "brackets" for 30 pseudoknot levels (and 19 more encoded with cyrillic letters)
    # Higher levels will be simply ignored
    levels = ['()','[]','{}','<>','Aa','Bb','Cc','Dd','Ee','Ff','Gg',
              'Hh','Ii','Jj','Kk','Ll','Mm','Nn','Oo','Pp','Qq','Rr',
              'Ss','Tt','Uu','Vv','Ww','Xx','Yy','Zz',
              'Бб','Гг','Дд','Ёё','Жж','Йй','Лл','Пп',
              'Фф','Цц','Чч','Шш','Щщ','Ьь','Ыы','Ъъ','Ээ','Юю','Яя']

    # groups of non-conflicting base pairs
    groups = [set(),]

    # normalize the pairs (i.e. ensure v < w)
    pairs = set((min(v, w), max(v, w)) for v, w in newpairs)
    
    for pair in sorted(pairs):

        level = 0

        # find the minimum level where the pair is not in conflict
        # with any base pair of that level
        while any(v[0]<=pair[0]<=v[1]<=pair[1] or
                  pair[0]<=v[0]<=pair[1]<=v[1] for v in groups[level]):
            level += 1
            if level == len(groups):
                groups.append(set())
            if level == len(levels):
                levels.append('..')

        # add the pair to the determined level
        groups[level].add(pair)

    # kind of a bubble sort of the base pairs among the levels
    # to maximize the number of base pairs of the lowest levels
    # e.g. to turn (..[[[...)...]]] into [..(((...]...)))
    for times in range(len(groups)-1):
        for i in range(len(groups)-1):

            rest = {v for v in groups[i+1] if any(v[0]<=w[0]<=v[1]<=w[1] or
                                                       w[0]<=v[0]<=w[1]<=v[1]
                                                       for w in groups[i])}
            clean = groups[i+1] - rest

            while rest:

                confjprev = set()
                confiprev = set()

                confj = rest.pop()
                rest.add(confj)
                confj = {confj,}
                confi = {v for v in groups[i] if any(v[0]<=w[0]<=v[1]<=w[1] or
                                                     w[0]<=v[0]<=w[1]<=v[1]
                                                     for w in confj)}

                while confjprev != confj or confiprev != confi:

                    confjprev = confj
                    confiprev = confi

                    confj = {v for v in rest if any(v[0]<=w[0]<=v[1]<=w[1] or
                                                    w[0]<=v[0]<=w[1]<=v[1]
                                                    for w in confi)}
                    confi = {v for v in groups[i] if any(v[0]<=w[0]<=v[1]<=w[1] or
                                                     w[0]<=v[0]<=w[1]<=v[1]
                                                     for w in confj)}

                if len(confi) < len(confj):

                    groups[i]   = confj | (groups[i] - confi)
                    groups[i+1] = confi | (groups[i+1] - confj)

                rest = rest - confj

            if clean:

                groups[i] |= clean
                groups[i+1] -= clean

    if returnlevels:
        levels = {}
        for lev, group in enumerate(groups):
            for bp in group:
                levels[bp] = lev + 1
        return levels

    # remove all levels higher than levellimit (if specified)
    if levellimit >= 0:
        groups = groups[:levellimit]

    # add all the pairs to the dbn string
    # according to their levels  
    for i, group in enumerate(groups):
        for pair in group:
            dbn[pair[0]] = levels[i][0]
            dbn[pair[1]] = levels[i][1]
            
    return ''.join(dbn)


def DBNToPairs(dbn):
    """Convert the dbn string into a sorted list of base pairs"""
    pairs = set()

    # keys == closing brackets, values == matching opening brackets
    closing = {'>':'<',']':'[',')':'(','}':'{','a':'A','b':'B','c':'C','d':'D',
               'e':'E','f':'F','g':'G','h':'H','i':'I','j':'J','k':'K','l':'L',
               'm':'M','n':'N','o':'O','p':'P','q':'Q','r':'R','s':'S','t':'T',
               'u':'U','v':'V','w':'W','x':'X','y':'Y','z':'Z',
               'б':'Б','г':'Г','д':'Д','ё':'Ё','ж':'Ж','й':'Й','л':'Л','п':'П',
               'ф':'Ф','ц':'Ц','ч':'Ч','ш':'Ш','щ':'Щ','ь':'Ь','ы':'Ы','ъ':'Ъ',
               'э':'Э','ю':'Ю','я':'Я',}
    # 30+19 bp stacks for 30+19 allowed pseudoknot levels
    stack = {'<':[],'(':[],'{':[],'[':[],'A':[],'B':[],'C':[],'D':[],'E':[],
             'F':[],'G':[],'H':[],'I':[],'J':[],'K':[],'L':[],'M':[],'N':[],
             'O':[],'P':[],'Q':[],'R':[],'S':[],'T':[],'U':[],'V':[],'W':[],
             'X':[],'Y':[],'Z':[],
             'Б':[],'Г':[],'Д':[],'Ё':[],'Ж':[],'Й':[],'Л':[],'П':[],'Ф':[],
             'Ц':[],'Ч':[],'Ш':[],'Щ':[],'Ь':[],'Ы':[],'Ъ':[],'Э':[],'Ю':[],
             'Я':[],}
              
    for i,v in enumerate(dbn):
        # if we observe an opening bracket
        # then add its index into the matching stack
        if v in stack: 
            stack[v].append(i)
        # else if we observe the closing bracket
        # take the opening index from the matching stack
        # and add the base pair to the pairs set
        elif v in closing:
            # this is to handle closing brackets with no
            # opening partner - they will be ignored
            if stack[closing[v]]:
                pairs.add((stack[closing[v]].pop(), i))

    return sorted(pairs)


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


def PredictSQUARNAs1(dataset, fam):

    command = "python SQUARNA.py i={} a step3=1 > outp3.tmp".format("datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictSQUARNAs2(dataset, fam):

    command = "python SQUARNA.py i={} a step3=2 > outp3.tmp".format("datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictSQUARNAs3i(dataset, fam):

    command = "python SQUARNA.py i={} a step3=i > outp3.tmp".format("datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictSQUARNAs3u(dataset, fam):

    command = "python SQUARNA.py i={} a step3=u > outp3.tmp".format("datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictRNAalifold(dataset, fam):

    command = "~/software/ViennaRNA-2.7.0/src/bin/RNAalifold --noPS {} > outp3.tmp".format("datasets/{}/aln/{}.aln".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictCentroidAlifold(dataset, fam):

    command = "~/software/centroid-rna-package-master/build/src/centroid_alifold"+\
              " {} > outp3.tmp".format("datasets/{}/aln/{}.aln".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictIPknot(dataset, fam):

    command = "~/software/ipknot-master/build/ipknot"+\
              " {} > outp3.tmp".format("datasets/{}/aln/{}.aln".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].strip()


def PredictRscapeNested(dataset, fam):

    os.makedirs("tmp",   exist_ok = True)
    command = "cd tmp; ~/software/rscape_v2.5.6/bin/R-scape"+\
              " --cacofold --covmin 4 --nofigures --rna --outname outp3.tmp"+\
              " {} > outp3.tmp".format("../datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    headers, seqnames, seqdict, gcnames, gcdict = ReadStockholm("tmp/outp3.tmp.cacofold.sto")
    return PairsToDBN(DBNToPairs(gcdict["SS_cons"]),len(gcdict["SS_cons"]))             
                 

def PredictRscapeTotal(dataset, fam):

    os.makedirs("tmp",   exist_ok = True)
    command = "cd tmp; ~/software/rscape_v2.5.6/bin/R-scape"+\
              " --cacofold --covmin 4 --nofigures --rna --outname outp3.tmp"+\
              " {} > outp3.tmp".format("../datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    headers, seqnames, seqdict, gcnames, gcdict = ReadStockholm("tmp/outp3.tmp.cacofold.sto")

    pairs = DBNToPairs(gcdict["SS_cons"])
    seen = set(pos for bp in pairs for pos in bp)

    for name in gcnames:
        if name.startswith("SS_cons_"):
            for v, w in DBNToPairs(gcdict[name]):
                if v not in seen and w not in seen:
                    seen.add(v)
                    seen.add(w)
                    pairs.append((v, w))
    return PairsToDBN(sorted(pairs), len(gcdict["SS_cons"])) 


if __name__ == "__main__":

    #dataset = "SubAli" # RNAStralignExt / Rfam14.9 / RfamPDB / SubAli / SeqSim
    #tool    = "IPknot"

    '''for dataset, tool in (("RNAStralignExt","SQUARNAs1"),
                          ("RNAStralignExt","SQUARNAs2"),
                          ("RNAStralignExt","SQUARNAs3i"),
                          ("RNAStralignExt","SQUARNAs3u"),
                          ("Rfam14.9","SQUARNAs1"),
                          ("Rfam14.9","SQUARNAs2"),
                          ("Rfam14.9","SQUARNAs3i"),
                          ("Rfam14.9","SQUARNAs3u"),
                          ("RfamPDB","SQUARNAs1"),
                          ("RfamPDB","SQUARNAs2"),
                          ("RfamPDB","SQUARNAs3i"),
                          ("RfamPDB","SQUARNAs3u"),
                          ("SubAli","SQUARNAs1"),
                          ("SubAli","SQUARNAs2"),
                          ("SubAli","SQUARNAs3i"),
                          ("SubAli","SQUARNAs3u"),
                          ("SeqSim","SQUARNAs1"),
                          ("SeqSim","SQUARNAs2"),
                          ("SeqSim","SQUARNAs3i"),
                          ("SeqSim","SQUARNAs3u"),
                          ("S01AliCMclean","SQUARNAs1"),
                          ("S01AliCMclean","SQUARNAs2"),
                          ("S01AliCMclean","SQUARNAs3i"),
                          ("S01AliCMclean","SQUARNAs3u"),
                          ("S01AliUngapclean","SQUARNAs1"),
                          ("S01AliUngapclean","SQUARNAs2"),
                          ("S01AliUngapclean","SQUARNAs3i"),
                          ("S01AliUngapclean","SQUARNAs3u"),
                          ("S01Aliclean","SQUARNAs1"),
                          ("S01Aliclean","SQUARNAs2"),
                          ("S01Aliclean","SQUARNAs3i"),
                          ("S01Aliclean","SQUARNAs3u"),):'''

    for dataset, tool in (("RNAStralignExt", "RscapeTotal"),
                          ):
                
        outname = "{}_{}".format(dataset,tool)
        title = '\t'.join("NAME LEN DEPTH TIME TP FP FN PRC RCL FS DBN PRED".split())
        outp1 = open(outname+'.fas','w')
        outp2 = open(outname+'.tsv','w')
        outp2.write(title+'\n')

        t0 = time.time()

        famfiles = glob.glob("datasets/{}/sto/*".format(dataset))
        
        fams = []

        for famfile in famfiles:
            fam = os.path.basename(famfile).split('.')[0]
            headers, seqnames, seqdict, gcnames, gcdict = ReadStockholm(famfile)
            LEN    = len(gcdict['SS_cons'])
            DEPTH  = len(seqnames)
            refdbn = PairsToDBN(DBNToPairs(gcdict['SS_cons']), LEN)
            fams.append((LEN, DEPTH, fam, refdbn))
        fams.sort()

        cnt = 0  
        for LEN, DEPTH, fam, refdbn in fams:
            cnt += 1
            name = '>'+fam
            print(name,end='')

            preddbn = { "SQUARNAs1":PredictSQUARNAs1,
                        "SQUARNAs2":PredictSQUARNAs2,
                        "SQUARNAs3i":PredictSQUARNAs3i,
                        "SQUARNAs3u":PredictSQUARNAs3u,
                        "RNAalifold":PredictRNAalifold,
                        "CentroidAlifold":PredictCentroidAlifold,
                        "IPknot": PredictIPknot,
                        "RscapeNested": PredictRscapeNested,
                        "RscapeTotal" : PredictRscapeTotal,
                        }[tool](dataset, fam)

            t1 = time.time()-t0

            print("...COMPLETE ({}sec) == {}/{}".format(round(t1,3), cnt, len(fams)))

            pairsr = set(DBNToPairs(refdbn))
            pairsq = set(DBNToPairs(preddbn))

            TP = len(pairsr & pairsq)
            FP = len(pairsq - pairsr)
            FN = len(pairsr - pairsq)
                    
            FS = 2*TP / (2*TP + FN + FP) if (TP + FN + FP) else 1
            PRC = (TP / (TP + FP)) if (TP+FP) else 1
            RCL = (TP / (TP + FN)) if (TP+FN) else 1

            outp1.write(name+'\n')
            outp1.write(refdbn+'\n')
            outp1.write(preddbn+'\n')
            outp1.write("LEN={} DEPTH={}, TIME={}sec TP={} FP={} FN={} PRC={} RCL={} FS={}\n"\
                        .format(LEN, DEPTH,round(t1,3),
                                TP,FP,FN,
                                round(PRC,3),
                                round(RCL,3),
                                round(FS,3)))
            res = [name[1:], LEN, DEPTH, round(t1,3), TP, FP, FN,
                    round(PRC,3), round(RCL,3), round(FS,3),
                    refdbn,preddbn]
            outp2.write('\t'.join([str(g) for g in res])+'\n')

        outp1.close()
        outp2.close()

