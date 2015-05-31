'''
    aPPRove: Preddicting PPR-RNA Interaction Using Primary Structure
    Copyright (C) 2015 Thomas Harrison

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
from multiprocessing import Process, Queue
import ppr_phmm,psScan
import argparse,os,sys
import json
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from scipy import stats

'''
This is the driver script of aPPRove. 
It splits the user input into jobs and runs them.
It will then output the results to stdout
'''


def getNumOverThreshold(scores):
    t=scores[-1]
    n=0.0
    m=0.0
    for i in range(len(scores)-1):
        m+=1
        if scores[i] >= t:
            n+=1    
    return n/m

def run(pairs,mPairs,outFile,dbScores):
    '''The Multithreaded part alligns aPPRove on all of the target transcripts and then
       Dumps each score into a json text file'''
    queueStr=""
    dumpList=[]
    for t in pairs:
        res=ppr_phmm.run(t[1],mPairs)
        name = args.proteinFile.split("/")[-1]
        queueStr="Bining Protein to Target : "+ t[0]+"\n"
        for i in range(len(res)):
            firstM=1
            for z in range (len(res[i][0][1])):
                if "m" in res[i][0][1][z]:
                    firstM=z+1
                    break
            queueStr+= "Alignment Starts at : "+str(firstM)+" on transcript. \n"
            queueStr+= "Alignment : "+str(i+1)+"\n"
            queueStr+= res[i][0][0]+"\n"
            queueStr+= "Has Score: "+str(res[i][1])+"\n"
            queueStr+= "Raw: " +str(res[i][2])+"\n"
            if not dbScores == None:
                dbScores.append(res[i][1])
                scores=stats.zscore(np.array(dbScores))
                fp = getNumOverThreshold(scores)
                p_values = stats.norm.sf(scores)
                queueStr+= "p-value: " +str(p_values[-1])+"\n"
                dbScores.pop(len(dbScores)-1)


        queueStr+="######################################################################"
        dumpList.append((t[0],queueStr))

    outf=open(outFile,"w")
    json.dump(dumpList,outf)
    outf.close()












if __name__ == "__main__":

    #-p input protein file
    #-m a motif file (optional)
    #-t input target file
    #-d database file (output from al2db)
    #--footprint if a binding footprint is used
    #-k number of optimal alignments to be returned
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-p', dest='proteinFile')
    parser.add_argument('-m', dest='motif')
    parser.add_argument('-t', dest='target')
    parser.add_argument('-k', nargs='?', dest='k')
    parser.add_argument('-d', nargs='?', dest='databaseFile')
    parser.add_argument('-n', nargs='?', dest='n')
    

    mPairs=None
    tString = None
    k=7
    n=4


    
    args = parser.parse_args()
    dbFile= args.databaseFile
    
    if(args.proteinFile == None and args.motif == None):
        print "Please provide an input protein sequence"
        sys.exit()
    if args.motif == None:
        mPairs=  psScan.runScan(args.proteinFile)
    else:
        mPairs = args.motif.split(",")


    if args.target == None:
        print "Please provide a RNA target sequence"
        sys.exit()


    if not args.k == None:
        k = int(args.k)

    if not args.n == None:
        n = int(args.n)

    # parameterize phmm from Barkan et al
    
    ppr_phmm.setParams(0.005,0.0440414508008,0.0000000001,0.5,k)
    
    # form target and name pairs to be run
    targets=[]
    names=[]
    if "." in args.target:
        inp = open(args.target,"r")
        for line in inp:
            if ">" in line:
                names.append(line[1:].rstrip())
                targets.append("")
            else:
                targets[-1]+=line.rstrip()
    else:
        targets.append(args.target.rstrip())
        names.append("User Specified Seq")

    pairs = zip(names,targets)

    count = 0
    segementedTargets=[]

    # assign a job to a process and run
    for i in range(n):
        segementedTargets.append([])
    for t in pairs:
        if count == n:
            count=0
        segementedTargets[count].append(t)
        count+=1


    processes = {}
    resultsCollection=[]
    count=0

    dbScores=None
    if not dbFile == None:
        dbFile= open(dbFile,"r")
        dbScores= json.load(dbFile)
    for t in segementedTargets:
        processes[count]=Process(target=run, args=(t,mPairs,"ap"+str(count)+".out",dbScores))
        count+=1


    for key in processes:
        processes[key].start()



    for key in processes:
        processes[key].join()


    dumpCatch=[]

    #ouput results
    for key in processes:
        infile = open("ap"+str(key)+".out","r")
        dumpCatch.append(json.load(infile))
        infile.close()
        os.remove("ap"+str(key)+".out")

    for l in dumpCatch:
        for instance in l:
            resultsCollection.append(instance)

    resultsCollection.sort()
    for item in resultsCollection:
        print item[1]

    








