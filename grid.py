from SQUARNA.SQRNdbnseq import SQRNdbnseq
import numpy as np
from multiprocessing import Pool
import os


def MainFunc(inputparams):

    queue,suboptmin,suboptmax,suboptsteps,bpp,minfinscorefactor,distcoef,orderpenalty,GC,AU,GU,minlen,minbpscore = inputparams
    
    paramsets = []
    paramsets.append({"bpweights" : {'GU' : GU,
                                    'AU' : AU,
                                    'GC' : GC,},
                        "suboptmax" : suboptmax,
                        "suboptmin" : suboptmin,
                        "suboptsteps": suboptsteps,
                        "minlen" : minlen,
                        "minbpscore" : minbpscore,
                        "minfinscorefactor" : minfinscorefactor,
                        "distcoef" : distcoef,
                        "bracketweight" :  -2,
                        "orderpenalty"  : orderpenalty,
                        "loopbonus": 0.125,
                        "maxstemnum" : maxstemnum,
                        "bpp": bpp,
                        })

    resultsC = []
    resultsB = []
    

    for obj in queue:

        name, seq, dbn, rst, react = obj

        result = SQRNdbnseq(seq, react, rst, dbn,
                            paramsets, conslim, toplim,
                            algos = algos,
                            threads = 1, mp = False)
                                                
        resultsC.append(result[2])
        resultsB.append(result[3])

    return inputparams, (resultsC, resultsB)
    

if __name__ == "__main__":

    from collections import Counter
  
    queue = []

    with open("datasets/SRtrain150.fas") as file:
        lines = file.readlines()

        for ii in range(0,len(lines)-2,3):

            nm  = lines[ii].strip()[1:]
            sq  = lines[ii+1].strip()
            db  = lines[ii+2].strip()
            queue.append([nm, sq, db, None, None])

    outpname = 'Greedy.tsv'

    outp  = open(outpname,'w')
    algos = set("G")

    title = '\t'.join("bpp minlen minbpscore GC AU GU minfinscorefactor distcoef orderpenalty tpc fpc fnc fstotc fsc prc rcc tp5 fp5 fn5 fstot5 fs5 pr5 rc5 harmon5".split())
    print(title)
    outp.write(title+'\n')
    outp.close()

    ######################################

    threads = os.cpu_count()

    conslim = 1
    toplim  = 5
    maxstemnum = 10**6

    inputs = []

    for suboptmin in (0.65,):
        for suboptmax in (0.9,):
            for suboptsteps in (1,):
                for bpp in (1, 0.5, 0, -0.5, -1):
                    for minfinscorefactor in (0.99, 1.25):
                        for distcoef in (0.09, 0.1):
                            for orderpenalty in (1.0, 1.35):
                                for GC in (2.0, 2.5, 3.0, 3.25, 3.5, 4.0):
                                    for AU in (0.5, 1.0, 1.25, 1.5, 2.0, 2.5):
                                        for GU in (-1.5, -1.25, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 1.5):
                                            for minlen,minbpscore in ((2, AU+GU),(2, AU+GU-0.25),(2, AU+GU+0.25),(2, AU+GU-0.5),(2, AU+GU+0.5),(2, AU+GU-0.75),(2, AU+GU+0.75),
                                                                      (2, AU+AU),(2, AU+AU-0.25),(2, AU+AU+0.25),(2, AU+AU-0.5),(2, AU+AU+0.5),(2, AU+AU-0.75),(2, AU+AU+0.75),
                                                                      (2, GC+GU),(2, GC+GU-0.25),(2, GC+GU+0.25),(2, GC+GU-0.5),(2, GC+GU+0.5),(2, GC+GU-0.75),(2, GC+GU+0.75),
                                                                      (2, GC+AU),(2, GC+AU-0.25),(2, GC+AU+0.25),(2, GC+AU-0.5),(2, GC+AU+0.5),(2, GC+AU-0.75),(2, GC+AU+0.75),
                                                                      (2, GC+GC),(2, GC+GC-0.25),(2, GC+GC+0.25),(2, GC+GC-0.5),(2, GC+GC+0.5),(2, GC+GC-0.75),(2, GC+GC+0.75),):

                                                if minbpscore < 0:
                                                    continue
                                                inputs.append([queue,suboptmin,suboptmax,suboptsteps,bpp,minfinscorefactor,distcoef,orderpenalty,GC,AU,GU,minlen,minbpscore])

    with Pool(threads) as pool:
        for inputparams,result in pool.imap(MainFunc, inputs):

            resultsC, resultsB = result
            queue,suboptmin,suboptmax,suboptsteps,bpp,minfinscorefactor,distcoef,orderpenalty,GC,AU,GU,minlen,minbpscore = inputparams                          
            print(bpp, minlen, minbpscore, GC, AU, GU, minfinscorefactor, distcoef, orderpenalty, sep='\t', end='\t')

            tpC = sum(x[0] for x in resultsC)
            fpC = sum(x[1] for x in resultsC)
            fnC = sum(x[2] for x in resultsC)
            fsC = [x[3] for x in resultsC]
            prC = [x[4] for x in resultsC]
            rcC = [x[5] for x in resultsC]

            print(tpC, fpC, fnC, round(2*tpC / (2*tpC + fpC + fnC), 3), round(np.mean(fsC), 3),
                    round(np.mean(prC), 3), round(np.mean(rcC), 3), sep='\t', end='\t')

            tpB = sum(x[0] for x in resultsB)
            fpB = sum(x[1] for x in resultsB)
            fnB = sum(x[2] for x in resultsB)
            fsB = [x[3] for x in resultsB]
            prB = [x[4] for x in resultsB]
            rcB = [x[5] for x in resultsB]
            rkB = [x[6] for x in resultsB]

            fstot5  = 2*tpB / (2*tpB + fpB + fnB)
            fsmean5 = np.mean(fsB) 
                                                                        
            print(tpB, fpB, fnB, round(fstot5, 3), round(fsmean5, 3),
                    round(np.mean(prB), 3), round(np.mean(rcB), 3), round(2*fsmean5*fstot5/(fstot5+fsmean5), 3), sep = '\t')

            outp = open(outpname,'a')
            toprint = '\t'.join([str(xx) for xx in [bpp, minlen, minbpscore, GC, AU, GU,
                                                    minfinscorefactor, distcoef, orderpenalty,
                                                    tpC, fpC, fnC,
                                                    round(2*tpC / (2*tpC + fpC + fnC), 3),
                                                    round(np.mean(fsC), 3),
                                                    round(np.mean(prC), 3),
                                                    round(np.mean(rcC), 3),
                                                    tpB, fpB, fnB,
                                                    round(fstot5, 3),
                                                    round(fsmean5, 3),
                                                    round(np.mean(prB), 3),
                                                    round(np.mean(rcB), 3),
                                                    round(2*fsmean5*fstot5/(fstot5+fsmean5), 3)]])
            outp.write(toprint+'\n')
            outp.close()
