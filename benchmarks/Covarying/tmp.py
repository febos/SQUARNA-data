import os, glob


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


def WriteStockholm(stkfile,headers,seqnames,seqdict,gcnames,gcdict):

    famname = os.path.basename(stkfile).split('_')[0]

    newseqnames = []

    for seqname in seqnames:

        newseqname = famname+'__'+seqname if not seqname.startswith(famname) else seqname
        newseqnames.append(newseqname)
        seqdict[newseqname] = seqdict[seqname]

    for i in range(len(headers)):
        if headers[i].startswith("#=GR "):
            grname = headers[i].strip().split()[1]
            if not grname.startswith(famname):
                headers[i] = headers[i].replace(grname,famname+'__'+grname)

    seqnames = newseqnames

    grlens = [len(' '.join(h.strip().split()[:-1])) for h in headers if h.startswith("#=GR ")]

    longest = max([4,] + [len(x) for x in seqnames] + [len(x)+len('#=GC ') for x in gcnames] + grlens) + 1

    if not headers or 'STOCKHOLM' not in headers[0]:
        headers.insert(0,'# STOCKHOLM 1.0\n')

    if headers[0] != '# STOCKHOLM 1.0\n':
        headers[0] = '# STOCKHOLM 1.0\n'

    with open(stkfile,'w') as file:

        for line in headers:
            if not line.startswith("#=GR "):
                file.write(line)

        file.write('\n')

        for name in seqnames:
            file.write(name + ' '*(longest-len(name)) + seqdict[name] + '\n')

        for line in headers:
            if line.startswith("#=GR "):
                file.write(line)
        
        for name in gcnames:
            file.write('#=GC ' + name + ' '*(longest-len(name)-len('#=GC ')) + gcdict[name] + '\n')

        file.write('//\n')



for tool in "CentroidAlifold RNAalifold RscapeNested RscapeTotal SQUARNAs3u".split():

    files = sorted(glob.glob("Rfam/*.sto"))
    
    print(tool)
    preds = {}

    with open("Rfam14.9_{}.tsv".format(tool)) as inp:
        for line in inp:
            if line.startswith("NAME"):
                continue
            arr = line.strip().split()
            preds[arr[0]] = arr[-1]

    for file in files:
        fam = file.split('/')[-1].split('.')[0]
        print(fam,end=' ')
        headers, seqnames, seqdict, gcnames, gcdict = ReadStockholm(file)
        outfile = file.replace("Rfam/",tool+'/')
        dbn = preds[fam]
        pairs = DBNToPairs(dbn)
        levels = PairsToDBN(pairs,len(dbn),returnlevels = True)
        mx = max(levels.values()) if levels else 1
        gcnames = ['SS_cons','RF',]
        if mx == 1:
            gcdict['SS_cons'] = preds[fam]
        else:
            for ii in range(mx):
                layer = PairsToDBN([pair for pair in pairs if levels[pair]==ii+1],len(dbn))
                if ii == 0:
                    gcdict['SS_cons'] = layer
                else:
                    token = 'SS_cons_{}'.format(ii)
                    gcnames.append(token)
                    gcdict[token] = layer
        WriteStockholm(outfile,headers, seqnames, seqdict, gcnames, gcdict)
    print()
        
