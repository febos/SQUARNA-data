import os, time

def CombinePairsToDBN(newpairs, length = 0, returnlevels = False, levellimit = -1):
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

def StemmedIsolated(sorted_pairs):

    sp = sorted_pairs

    stemmed, isolated = [], []

    if len(sp) < 2:
        if not sp:
            return [], []
        else:
            return [],[sp[0][0],sp[0][1]]

    for i in range(len(sp)):

        if i==0:
            if sp[i][0] + 1 == sp[i+1][0] and sp[i][1] == sp[i+1][1] + 1:
                stemmed.append(sp[i])
            else:
                isolated.append(sp[i][0])
                isolated.append(sp[i][1])
                
        elif i==len(sp)-1:
            if sp[i-1][0] + 1 == sp[i][0] and sp[i-1][1] == sp[i][1] + 1:
                stemmed.append(sp[i])
            else:
                isolated.append(sp[i][0])
                isolated.append(sp[i][1])
            
        else:
            if sp[i][0] + 1 == sp[i+1][0] and sp[i][1] == sp[i+1][1] + 1 or \
               sp[i-1][0] + 1 == sp[i][0] and sp[i-1][1] == sp[i][1] + 1:
                stemmed.append(sp[i])
            else:
                isolated.append(sp[i][0])
                isolated.append(sp[i][1])

    return stemmed, isolated  

def GetPairs(dbn):

    pairs = set()

    closing = {'>':'<',']':'[',')':'(','}':'{','a':'A','b':'B','c':'C','d':'D','e':'E'}
    stack = {'<':[],'(':[],'{':[],'[':[],'A':[],'B':[],'C':[],'D':[],'E':[]}

    for i,v in enumerate(dbn):

        if v in stack:
            stack[v].append(i)
        if v in closing:
            if stack[closing[v]]:
                pairs.add((stack[closing[v]].pop(),i))

    return sorted(pairs)
          

def NoLone(dbn):

    pairs = GetPairs(dbn)
    lone_pos_set = set(StemmedIsolated(pairs)[1])

    return ''.join([dbn[i] if i not in lone_pos_set else '.' for i in range(len(dbn))])


def PredictRNAfold(seq):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')

    os.system("~/software/ViennaRNA-2.7.0/src/bin/RNAfold --noPS < inp.tmp > outp2.tmp")
    with open("outp2.tmp") as outp:
        dbn = outp.readlines()[2].split()[0]
    return [dbn,]

def PredictRNAsubopt5(seq, top = 5):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')

    if len(seq) > 2000:
        add = "--deltaEnergy=0.1"
    else:
        add = ""
    os.system("~/software/ViennaRNA-2.7.0/src/bin/RNAsubopt --sorted {} < inp.tmp > outp2.tmp".format(add))
    with open("outp2.tmp") as outp:
        dbns = [x.split()[0] for x in outp.readlines()[2:]]
    return dbns[:top]


def PredictIPknot(seq):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')

    os.system("~/software/ipknot-master/build/ipknot inp.tmp > outp2.tmp")
    with open("outp2.tmp") as outp:
        dbn = outp.readlines()[2].strip()
    return [dbn,] 

def PredictMXfold2(seq):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')

    os.system("mxfold2 predict inp.tmp > outp2.tmp")
    with open("outp2.tmp") as outp:
        dbn = outp.readlines()[2].split()[0]
    return [dbn,] 

def BPSEQtoDBN(bpseqfile):

    with open(bpseqfile) as outp:
        pairs = []
        lines = outp.readlines()[1:]
        for line in lines:
            x,y,z = line.strip().split()
            x = int(x)
            z = int(z)
            if x < z:
                pairs.append((x-1,z-1))
        return CombinePairsToDBN(pairs,x)


def PredictCONTRAfold(seq):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')

    os.system("~/software/contrafold/src/contrafold predict inp.tmp > outp2.tmp")
    with open("outp2.tmp") as outp:
        dbn = outp.readlines()[3].strip()
    return [dbn,] 

def PredictEternaFold(seq):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')

    params = "--params ~/software/EternaFold-1.3.1/parameters/EternaFoldParams.v1"
    os.system("~/software/contrafold/src/contrafold predict inp.tmp {} > outp2.tmp".format(params))
    with open("outp2.tmp") as outp:
        dbn = outp.readlines()[3].strip()
    return [dbn,] 


def PredictSPOTRNA(seq):

    seqs = [seq,]

    inpf = "inp.tmp"
    names = []
    with open(inpf,'w') as inp:
        for ii,seq in enumerate(seqs):
            names.append("seq%s"%ii)
            inp.write('>seq%s\n'%ii)
            inp.write(seq+'\n')

    os.system('conda run -n spotrna python ~/software/SPOT-RNA/SPOT-RNA.py --inputs inp.tmp --outputs "tmp" --cpu 32')

    res = []

    for name in names:
        res.append(BPSEQtoDBN("tmp/{}.bpseq".format(name)))
    return res

def CTtoDBN(ct_file):

    res = []

    with open(ct_file) as file:
        lines = file.readlines()

    while lines:

        ln = int(lines[0].split()[0])

        mfe = lines[1:ln+1]

        skpairs = set()
        
        for line in mfe:
            linesplit = line.split()
            pair = (int(linesplit[0]) - 1,int(linesplit[4]) - 1)

            if pair[-1] == -1 or not pair[0] < pair[1]:
                continue
            skpairs.add(pair)

        res.append(CombinePairsToDBN(sorted(skpairs),len(mfe)))
        lines = lines[ln+1:]

    return res


def PredictShapeKnots(seq, top = 1):

    SHAPEKNOTS_PATH   = "~/software/RNAstructure6.5/exe/ShapeKnots-smp"
    SHAPEKNOTS_TABLES = "DATAPATH=~/software/RNAstructure6.5/data_tables"

    with open("inp.tmp","w") as inp:
        inp.write(">seq"+'\n')
        inp.write(seq+'\n')

    os.system("{} {} inp.tmp outp2.tmp".format(SHAPEKNOTS_TABLES, SHAPEKNOTS_PATH))

    try:
        res = CTtoDBN("outp2.tmp")[:top]
        if not res:
            return ['.'*len(seq),]
        return res
    except:
        return ['.'*len(seq),]


def PredictShapeKnots5(seq):
    return PredictShapeKnots(seq, top = 5)


def PredictSQUARNA(seq, conf = "def.conf", top = 1):

    with open("inp.tmp","w") as inp:
        inp.write(">seq"+'\n')
        inp.write(seq+'\n')

    if len(seq) >= 500:
        conf = "500.conf"
    if len(seq) >= 1000:
        conf = "1000.conf"
    print(', '+conf)
    
    os.system("python3 SQUARNA/SQUARNA.py i=inp.tmp c={} toplim={} > outp2.tmp".format(conf, top))

    cnt = 0
    flag = False
    res = []

    with open("outp2.tmp") as outp:
        for line in outp:

            if line.startswith("="):
                flag = True
                continue

            if flag:
                res.append(line.strip().split()[0])

    return res[:top]


def PredictSQUARNAE(seq):
    return PredictSQUARNA(seq, conf="edmonds.conf")

def PredictSQUARNAG(seq):
    return PredictSQUARNA(seq, conf="greedy.conf")

def PredictSQUARNAG5(seq):
    return PredictSQUARNA(seq, conf="greedy.conf", top = 5)

def PredictSQUARNAE5(seq):
    return PredictSQUARNA(seq, conf="edmonds.conf", top = 5)

def PredictSQUARNAGE(seq):
    return PredictSQUARNA(seq, algos="ge")

def PredictSQUARNAGE5(seq):
    return PredictSQUARNA(seq, top = 5, algos="ge")

def PredictSQUARNAn(seq):
    return PredictSQUARNA(seq, conf="nussinov.conf")

def PredictSQUARNAn5(seq):
    return PredictSQUARNA(seq, top = 5, algos="n")

def PredictSQUARNAGn(seq):
    return PredictSQUARNA(seq, algos="gn")

def PredictSQUARNAGn5(seq):
    return PredictSQUARNA(seq, top = 5, algos="gn")

def PredictSQUARNAghn(seq):
    return PredictSQUARNA(seq, algos="ghn")

def PredictSQUARNAghn5(seq):
    return PredictSQUARNA(seq, top = 5, algos="ghn")

def PredictSQUARNAghe(seq):
    return PredictSQUARNA(seq, algos="ghe")

def PredictSQUARNAghe5(seq):
    return PredictSQUARNA(seq, top = 5, algos="ghe")

def PredictSQUARNAgne(seq):
    return PredictSQUARNA(seq, algos="gne")

def PredictSQUARNAgne5(seq):
    return PredictSQUARNA(seq, top = 5, algos="gne")

def PredictSQUARNAghne(seq):
    return PredictSQUARNA(seq, algos="ghne")

def PredictSQUARNAghne5(seq):
    return PredictSQUARNA(seq, top = 5, algos="ghne")

def PredictSQUARNAhne(seq):
    return PredictSQUARNA(seq, algos="hne")

def PredictSQUARNAhne5(seq):
    return PredictSQUARNA(seq, top = 5, algos="hne")

def PredictSQUARNAh(seq):
    return PredictSQUARNA(seq, conf="hungarian.conf")

def PredictSQUARNAh5(seq):
    return PredictSQUARNA(seq, top = 5, algos="h")

def PredictSQUARNAGh(seq):
    return PredictSQUARNA(seq, algos="gh")

def PredictSQUARNAGh5(seq):
    return PredictSQUARNA(seq, top = 5, algos="gh")

def PredictSQUARNA5(seq):
    return PredictSQUARNA(seq, top = 5)

def PredictSQUARNAN(seq):
    return PredictSQUARNA(seq, top = 10**6)

def PredictSQUARNAalt(seq):
    return PredictSQUARNA(seq, conf = "alt.conf")

def PredictSQUARNAalt5(seq):
    return PredictSQUARNA(seq, conf = "alt.conf", top = 5)

def PredictSQUARNAnew(seq):
    return PredictSQUARNA(seq, conf = "newdef.conf")

def PredictSQUARNAnew5(seq):
    return PredictSQUARNA(seq, conf = "newdef.conf", top = 5)

def PredictSQUARNAaltN(seq):
    return PredictSQUARNA(seq, conf = "alt.conf", top = 10**6)

def PredictSQUARNAsk(seq):
    return PredictSQUARNA(seq, conf = "sk.conf")

def PredictSQUARNAsk5(seq):
    return PredictSQUARNA(seq, conf = "sk.conf", top = 5)

def PredictSQUARNAskN(seq):
    return PredictSQUARNA(seq, conf = "sk.conf", top = 10**6)

     

if __name__ == "__main__":

    NL      =  False
    
    dtst  = "SRtrain150"
    tl    = "CONTRAfold"

    for dataset, tool in (("SRtrain150", "SQUARNA"),
                          ):

        if NL:
            dataset += "NL"

        with open('datasets/{}.fas'.format(dataset)) as file:

            outname = "{}_{}".format(dataset,tool)
            title = '\t'.join("NAME LEN TIME RANK TP FP FN PRC RCL FS SEQ DBN PRED".split())
            outp1 = open(outname+'.fas','w')
            outp2 = open(outname+'.tsv','w')
            outp2.write(title+'\n')
            lines = file.readlines()

            t0 = time.time()
            
            for i in range(0,len(lines)-2,3):

                name = lines[i].strip()
                print(name,end='')
                seq = lines[i+1].strip().upper()
                dbn = lines[i+2].strip()

                structs = {"RNAfold":PredictRNAfold,
                           "IPknot": PredictIPknot,
                           "MXfold2":PredictMXfold2,
                           "SPOT-RNA":PredictSPOTRNA,
                           "SQUARNA": PredictSQUARNA,
                           "SQUARNA5": PredictSQUARNA5,
                           "SQUARNAnew": PredictSQUARNAnew,
                           "SQUARNAnew5": PredictSQUARNAnew5,
                           "SQUARNAE": PredictSQUARNAE,
                           "SQUARNAG": PredictSQUARNAG,
                           "SQUARNAG5": PredictSQUARNAG5,
                           "SQUARNAE5": PredictSQUARNAE5,
                           "SQUARNAGE": PredictSQUARNAGE,
                           "SQUARNAGE5": PredictSQUARNAGE5,
                           "SQUARNAn": PredictSQUARNAn,
                           "SQUARNAn5": PredictSQUARNAn5,
                           "SQUARNAGn": PredictSQUARNAGn,
                           "SQUARNAGn5": PredictSQUARNAGn5,
                           "SQUARNAh": PredictSQUARNAh,
                           "SQUARNAh5": PredictSQUARNAh5,
                           "SQUARNAGh": PredictSQUARNAGh,
                           "SQUARNAGh5": PredictSQUARNAGh5,
                           "SQUARNAghn": PredictSQUARNAghn,
                           "SQUARNAghn5": PredictSQUARNAghn5,
                           "SQUARNAghe": PredictSQUARNAghe,
                           "SQUARNAghe5": PredictSQUARNAghe5,
                           "SQUARNAgne": PredictSQUARNAgne,
                           "SQUARNAgne5": PredictSQUARNAgne5,
                           "SQUARNAghne": PredictSQUARNAghne,
                           "SQUARNAghne5": PredictSQUARNAghne5,
                           "SQUARNAhne": PredictSQUARNAhne,
                           "SQUARNAhne5": PredictSQUARNAhne5,
                           "SQUARNAN": PredictSQUARNAN,
                           "SQUARNAalt": PredictSQUARNAalt,
                           "SQUARNAalt5": PredictSQUARNAalt5,
                           "SQUARNAaltN": PredictSQUARNAaltN,
                           "SQUARNAsk": PredictSQUARNAsk,
                           "SQUARNAsk5": PredictSQUARNAsk5,
                           "SQUARNAskN": PredictSQUARNAskN,
                           "ShapeKnots": PredictShapeKnots,
                           "ShapeKnots5": PredictShapeKnots5,
                           "RNAsubopt5": PredictRNAsubopt5,
                           "CONTRAfold": PredictCONTRAfold,
                           "EternaFold": PredictEternaFold,
                           }[tool](seq)

                t1 = time.time()-t0

                print("...COMPLETE ({}sec)".format(round(t1,3)))

                

                # Clean <3nt hairpins and non-canonical pairs
                structs = [CombinePairsToDBN([(v,w) for v,w in GetPairs(_)
                                              if w-v >= 4 and seq[v]+seq[w] in {'GC','CG',
                                                                                'GU','UG',
                                                                                'AU','UA'}],
                                             len(seq))
                           for _ in structs]

                # Clean lone bps if NL
                if NL:
                    structs = [NoLone(_) for _ in structs]

                best_ind    = -1
                best_fscore = -1
                BTP, BFP, BFN = -1, -1, -1

                pairsr = set(GetPairs(dbn))
                
                for i,pred in enumerate(structs):
                
                    pairsq = set(GetPairs(pred))

                    TP = len(pairsr & pairsq)
                    FP = len(pairsq - pairsr)
                    FN = len(pairsr - pairsq)
                
                    FS = 2*TP / (2*TP + FN + FP) if (TP + FN + FP) else 1
                    PRC = (TP / (TP + FP)) if (TP+FP) else 1
                    RCL = (TP / (TP + FN)) if (TP+FN) else 1

                    if FS > best_fscore:
                        best_ind = i
                        best_fscore = FS
                        BTP, BFP, BFN = TP, FP, FN

                best_prc = (BTP / (BTP + BFP)) if (BTP+BFP) else 1
                best_rcl = (BTP / (BTP + BFN)) if (BTP+BFN) else 1

                outp1.write(name+'\n')
                outp1.write(seq+'\n')
                outp1.write(dbn+'\n')
                for pred in structs:
                    outp1.write(pred+'\n')
                outp1.write("LEN={} TIME={}sec RANK={} TP={} FP={} FN={} PRC={} RCL={} FS={}\n"\
                            .format(len(seq),round(t1,3),best_ind+1,
                                    BTP,BFP,BFN,
                                    round(best_prc,3),
                                    round(best_rcl,3),
                                    round(best_fscore,3)))
                
                res = [name[1:], len(seq), round(t1,3), best_ind+1, BTP,BFP,BFN,
                       round(best_prc,3), round(best_rcl,3), round(best_fscore,3),
                       seq,dbn,structs[best_ind]]
                outp2.write('\t'.join([str(g) for g in res])+'\n')

        outp1.close()
        outp2.close()


