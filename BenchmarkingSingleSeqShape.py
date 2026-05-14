import os, time, glob

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


def PredictRNAfold(seq, react):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')
    with open("inp.dat",'w') as inp:
        for k,rval in enumerate(react.split()):
            inp.write("{} {}\n".format(k+1,rval))

    os.system("~/software/ViennaRNA-2.7.1/src/bin/RNAfold --shape=inp.dat --noPS < inp.tmp > outp2.tmp")
    with open("outp2.tmp") as outp:
        dbn = outp.readlines()[2].split()[0]
    return [dbn,]       


def PredictRNAsubopt5(seq, react, top = 5):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')
    with open("inp.dat",'w') as inp:
        for k,rval in enumerate(react.split()):
            inp.write("{} {}\n".format(k+1,rval))

    os.system("~/software/ViennaRNA-2.7.1/src/bin/RNAsubopt --shape=inp.dat --sorted < inp.tmp > outp2.tmp")
    with open("outp2.tmp") as outp:
        dbns = [x.split()[0] for x in outp.readlines()[2:]]
    return dbns[:top]     


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


def PredictShapeKnots(seq, react, top = 1):

    SHAPEKNOTS_PATH   = "~/software/RNAstructure6.5/exe/ShapeKnots-smp"
    SHAPEKNOTS_TABLES = "DATAPATH=~/software/RNAstructure6.5/data_tables"

    with open("inp.tmp","w") as inp:
        inp.write(">seq"+'\n')
        inp.write(seq+'\n')
    with open("inp.dat",'w') as inp:
        for k,rval in enumerate(react.split()):
            inp.write("{} {}\n".format(k+1,rval))

    os.system("{} {} inp.tmp outp2.tmp -sh inp.dat".format(SHAPEKNOTS_TABLES,
                                                           SHAPEKNOTS_PATH))

    try:
        res = CTtoDBN("outp2.tmp")[:top]
        if not res:
            return ['.'*len(seq),]
        return res
    except:
        print('FAILED')
        return ['.'*len(seq),]


def PredictShapeKnots5(seq, react):
    return PredictShapeKnots(seq, react, top = 5)


def PredictShapify(seq, react):

    #inpf = "inp.tmp"
    #with open(inpf,'w') as inp:
    #    inp.write('>seq\n')
    #    inp.write(seq+'\n')
    with open("inp.dat",'w') as inp:
        for rval in react.split():
            inp.write("{}\n".format(rval))

    os.system('~/software/Shapify/HFold_iterative --s "{}" --shape "inp.dat" > outp2.tmp'.format(seq))
    with open("outp2.tmp") as outp:
        lines = outp.readlines()
        try:
            theline = [x for x in lines if x.startswith('Result_0:')][0]
        except Exception as err:
            print(lines)
            raise err
        dbn = theline.strip().split()[1]
    return [dbn,]

def PredictShapify5(seq, react):

    #inpf = "inp.tmp"
    #with open(inpf,'w') as inp:
    #    inp.write('>seq\n')
    #    inp.write(seq+'\n')
    with open("inp.dat",'w') as inp:
        for rval in react.split():
            inp.write("{}\n".format(rval))

    os.system('~/software/Shapify/HFold_iterative --s "{}" --shape "inp.dat" --n 5 > outp2.tmp'.format(seq))
    with open("outp2.tmp") as outp:
        lines = outp.readlines()
        thelines = [x for x in lines if x.startswith('Result')]
        if thelines:
            dbns = [theline.strip().split()[1] for theline in thelines][:5]
        else:
            dbns = [lines[1].strip().split()[0],]
        
    return dbns

def PredictSQUARNA(seq, react, conf = "def.conf", top = 1):

    with open("inp.tmp","w") as inp:
        inp.write(">seq"+'\n')
        inp.write(seq+'\n')
        inp.write(react+'\n')

    if len(seq) >= 500:
        conf = "500.conf"
    if len(seq) >= 1000:
        conf = "1000.conf"
    print(', '+conf)
    if conf in ('def.conf','500.conf','1000.conf'):
        os.system("SQUARNA i=inp.tmp toplim={} > outp2.tmp".format(top))
    else:
        os.system("SQUARNA i=inp.tmp c={} toplim={} > outp2.tmp".format(conf,top))

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


def PredictSQUARNA5(seq, react):
    return PredictSQUARNA(seq, react, top = 5)

def PredictSQUARNAN(seq, react):
    return PredictSQUARNA(seq, react, conf = "nussinov.conf")

def PredictSQUARNAalt(seq, react, reactfile):
    return PredictSQUARNA(seq, react, reactfile, conf = "alt.conf")

def PredictSQUARNAalt5(seq, react, reactfile):
    return PredictSQUARNA(seq, react, reactfile, conf = "alt.conf", top = 5)

def PredictSQUARNAaltN(seq, react, reactfile):
    return PredictSQUARNA(seq, react, reactfile, conf = "alt.conf", top = 10**6)

def PredictSQUARNAsk(seq, react, reactfile):
    return PredictSQUARNA(seq, react, reactfile, conf = "sk.conf")

def PredictSQUARNAsk5(seq, react, reactfile):
    return PredictSQUARNA(seq, react, reactfile, conf = "sk.conf", top = 5)

def PredictSQUARNAskN(seq, react, reactfile):
    return PredictSQUARNA(seq, react, reactfile, conf = "sk.conf", top = 10**6)

     

if __name__ == "__main__":
    
    dtst  = "S01"
    tl    = "ShapeKnots"

    for dataset, tool in (("S01",  "RNAfold"),
                          ("S01",  "SQUARNA"),
                          ):

        with open('datasets/{}.fas'.format(dataset)) as file:

            outname = "{}_{}".format(dataset,tool)
            title = '\t'.join("NAME LEN TIME RANK TP FP FN PRC RCL FS SEQ DBN PRED".split())
            outp1 = open(outname+'.fas','w')
            outp2 = open(outname+'.tsv','w')
            outp2.write(title+'\n')
            lines = file.readlines()

            t0 = time.time()
            
            for i in range(0,len(lines)-3,4):

                name = lines[i].strip()
                print(name,end='')
                seq = lines[i+1].strip().upper()
                dbn = lines[i+2].strip()
                react = lines[i+3].strip()
                #print(react)

                structs = {"RNAfold":PredictRNAfold,
                           "ShapeKnots": PredictShapeKnots,
                           "ShapeKnots5": PredictShapeKnots5,
                           "RNAsubopt5": PredictRNAsubopt5,
                           "Shapify": PredictShapify,
                           "Shapify5": PredictShapify5,
                           "SQUARNA": PredictSQUARNA,
                           "SQUARNA5": PredictSQUARNA5,
                           "SQUARNAN": PredictSQUARNAN,
                           "SQUARNAalt": PredictSQUARNAalt,
                           "SQUARNAalt5": PredictSQUARNAalt5,
                           "SQUARNAaltN": PredictSQUARNAaltN,
                           "SQUARNAsk": PredictSQUARNAsk,
                           "SQUARNAsk5": PredictSQUARNAsk5,
                           "SQUARNAskN": PredictSQUARNAskN,
                           }[tool](seq, react)

                t1 = time.time()-t0

                print("...COMPLETE ({}sec)".format(round(t1,3)))

                

                # Clean <3nt hairpins and non-canonical pairs
                structs = [CombinePairsToDBN([(v,w) for v,w in GetPairs(_)
                                              if w-v >= 4 and (seq[v]+seq[w]).upper().replace('T','U') in {'GC','CG',
                                                                                'GU','UG',
                                                                                'AU','UA'}],
                                             len(seq))
                           for _ in structs]

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


