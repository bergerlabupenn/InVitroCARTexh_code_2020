#!/usr/bin/python
# Charly noticed that there is an odd set of genes in her volcano plots which
# sit on a parabola and upon examination in the normalized DESeq2 counts, she
# sees that many of these have zero across the board except in one sample.
# This script filters out the genes in question using crude mathematical tools
# Greg Donahue, 07-30-2019

import sys, math

#counts_file = "Adjusted Counts (Only Positives).txt"
#deseq_file = "CS vs Day 0 Results.txt"
logfc_threshold = (-1.0,1.0)
logp_threshold = -1*math.log(0.05, 10)

usage = "USAGE: python find_extreme_genes.py ADJUSTED_COUNTS_TABLE "+ \
        "DESEQ_RESULTS"

def main(args):

    # Read command-line arguments
    if len(args) != 3: sys.exit(usage)
    counts_file = args[1]
    deseq_file = args[2]

    # Loads DESeq2 adjusted counts (function counts(dds, normalized=TRUE) in R)
    data = dict()
    f = open(counts_file); lines = f.readlines(); f.close()
    headers = lines[0].replace("\r", "")
    for line in lines[1:]:
        t = line.replace("\n", "").replace("\r", "").split("\t")
        data[t[0]] = [ float(x) for x in t[1:] ]

    # Loads DESeq2 results
    deseq = dict()
    f = open(deseq_file); lines = f.readlines(); f.close()
    for line in lines[1:]:
        t = line.replace("\n", "").replace("\r", "").split("\t")
        try: deseq[t[0]] = (float(t[2]),-1*math.log(float(t[-1]), 10))
        except Exception as e: deseq[t[0]] = (0.0,0.0)

    # Pulls out a gene list for genes up- or down-regulated
    results_up, results_down = list(), list()
    for k in list(deseq.keys()):
        if deseq[k][0] < logfc_threshold[0] and deseq[k][1] > logp_threshold:
            results_down.append(k)
        elif deseq[k][0] > logfc_threshold[1] and deseq[k][1] > logp_threshold:
            results_up.append(k)
            
    # Determined from IL22 and WDR63, negative arm
    c2 = (-15.7092,-1*math.log(0.000294742,10))
    c1 = (-22.69076461,-1*math.log(2.25219650758296E-08,10))
    slope = (c1[1]-c2[1])/(c1[0]-c2[0])
    intercept = c1[1]-c1[0]*slope
    bogus_down = list()
    print(headers)
    for k in results_down:
        if deseq[k][0] >= -10: continue
        D = getSupportVectorDistance(deseq[k], slope, intercept)
        if D < 1: bogus_down.append(k)
    for B in bogus_down: print(B+"\t"+"\t".join([ str(x) for x in data[B] ]))
    # Determined from ALK and INHBE, positive arm
    c1 = (22.3741316,-1*math.log(4.18923382767297E-08,10))
    c2 = (18.25474442,-1*math.log(0.0000173703118056835,10))
    slope = (c1[1]-c2[1])/(c1[0]-c2[0])
    intercept = c1[1]-c1[0]*slope
    bogus_up = list()
    for k in results_up:
        if deseq[k][0] <= 10: continue
        D = getSupportVectorDistance(deseq[k], slope, intercept)
        if D < 1: bogus_up.append(k)
    for B in bogus_up: print(B+"\t"+"\t".join([ str(x) for x in data[B] ]))

    # Check all genes not on the smile
    histogram = dict()
    for k in list(set(list(data.keys())).difference(set(bogus_down+bogus_up))):
        try: z = deseq[k]
        except Exception as e: continue
        if (deseq[k][0] < 1 and deseq[k][0] > -1) or deseq[k][1] < -1*math.log(0.05, 10): continue 
        try: histogram[countZeroes(data[k])] += 1
        except Exception as e: histogram[countZeroes(data[k])] = 1
    for k in sorted(list(histogram.keys())): print(str(k)+"\t"+str(histogram[k]))
    
    # Rewrite DESeq table with just filtered in genes
    f = open(deseq_file); lines = f.readlines(); f.close()
    f = open(deseq_file[:-3]+"Corrected.txt", 'w')
    g = open(deseq_file[:-3]+"Corrected_Total.txt", 'w')
    f.write(lines[0])
    g.write(lines[0])
    for line in lines[1:]:
        t = line[:-1].split("\t")
        if len(set([ t[0] ]).intersection(set(bogus_down+bogus_up))) > 0: continue
        g.write(line)
        try: z = deseq[t[0]]
        except Exception as e: continue
        if (deseq[t[0]][0] < 1 and deseq[t[0]][0] > -1) or deseq[t[0]][1] < -1*math.log(0.05, 10): continue
        f.write(line)
    f.close()

def countZeroes(vector):
    ret = 0
    for V in vector[:16]:
        if V == 0: ret += 1
    return ret

def getSupportVectorDistance(point, slope, intercept):
    y_low = point[0]*slope+intercept
    x_high = (point[1]-intercept)/slope
    T1 = math.fabs(point[1]-y_low)
    T2 = math.fabs(x_high-point[0])
    if T2 == 0 or T1 == 0: return 0.0
    A = math.atan(T1/T2)
    return math.sin(A)*T2

if __name__ == "__main__": main(sys.argv)
