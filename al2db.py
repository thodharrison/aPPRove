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

import matplotlib
# so it can be run as a background process
matplotlib.use('Agg')
from multiprocessing import Process, Queue
import ppr_phmm,psScan
import argparse,os,sys
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

'''
al2db builds up a distribution that can be fed into approve. 
This would allow for a user to get a statisticle context of an alignment to their own transcript. 
'''


def run(pairs,mPairs,outFile,s, q):

    '''The Multithreaded part alligns aPPRove on all of the chopped up transcripts and then
       Dumps each score into a json text file''' 
    queueStr=""
    dumpList=[]
    for t in pairs:
        res=ppr_phmm.run(t[1],mPairs)
        name=t[0]
        startIndex=int(t[0].split("@")[1])
        mStart=res[0][0][1].find("m")
        startIndex+=mStart

        normalizedScore=res[0][1]
        if (not betterFile == None and normalizedScore >= s):
            # We dont want to clog the pipe wand we dont want the file to be to big for the webserver
            if q.qsize() > 1500:
                q.put("Warning: Queue for better allignemnts is full. No more can be added.")

            else:
                q.put(t[0]+"\n")
                q.put(res[0][0][0]+"\n")
                q.put(str(normalizedScore)+"\n")
                print q.qsize()
            
        queueStr=t[0].split("@")[0]+"@"+str(startIndex)+"$"+str(normalizedScore)
        #print queueStr
        dumpList.append((t[0],queueStr))

    outf=open(outFile,"w")
    json.dump(dumpList,outf)
    outf.close()






if __name__ == "__main__":
    #-p protein file
    #-d database file
    #-o where output should be directed
    #-n no cores used
    #-s an optional score cutoff
    plt.ioff()
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-p', dest='proteinFile')
    parser.add_argument('-d', dest='target')
    parser.add_argument('-o', dest="outfiles")
    parser.add_argument('-n', nargs='?', dest='n')
    parser.add_argument('-s', nargs='?', dest='s')

    n=4
    s=0
    args = parser.parse_args()
    if(args.proteinFile == None):
        print "Please provide an input protein sequence"
        sys.exit()

    mPairs=  psScan.runScan(args.proteinFile)

    if args.target == None:
        print "Please provide a RNA target sequence"
        sys.exit()

    if not args.n == None:
        n = int(args.n)

    betterFile=None
    if not args.s == None:
        s = float(args.s)
        betterFile = open(args.outfiles+"_higherScores.txt","w")


    name=[]
    transcripts=[]

    inp = open(args.target,"r")
    for t in inp:
        if ">" in t:
            name.append(t[1:].rstrip())
            transcripts.append("")
        else:
            transcripts[-1]+=t.rstrip()

    # chop up the transcripts
    transcriptsChoppedNames=[]
    transcriptsChoppedSeqs=[]
    q=Queue()
    for z in range(len(name)):
        for i in range(len(transcripts[z])- len(mPairs)):
            transcriptsChoppedNames.append(name[z]+"@"+str(i))
            transcriptsChoppedSeqs.append(transcripts[z][i:i+len(mPairs)])



    #Parameterize phmm based off Barkan et al.

    ppr_phmm.setParams(0.005,0.0440414508008,0.0000000001,0.5,1)
    

    #zip names and sequences
    pairs = zip(transcriptsChoppedNames,transcriptsChoppedSeqs)
    count = 0
    segementedTargets=[]

    # start forming lists for each precessor
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
    # start all processes
    for t in segementedTargets:
        processes[count]=Process(target=run, args=(t,mPairs,str(count)+".out",s,q))
        count+=1


    for key in processes:
        processes[key].start()



    for key in processes:
        processes[key].join()


    dumpCatch=[]

    #ouput results and clean up
    for key in processes:
        infile = open(str(key)+".out","r")
        dumpCatch.append(json.load(infile))
        infile.close()
        os.remove(str(key)+".out")

    for l in dumpCatch:
        for instance in l:
            resultsCollection.append(instance)

    resultsCollection.sort()

    uniques={}
    for item in resultsCollection:
        #print item[1]
        uniques[item[1].split("$")[0]] =float( item[1].split("$")[1])

    forHist=[]
    for item in uniques:
         forHist.append(uniques[item])

    plt.hist(forHist, bins=len(set(forHist)))
    plt.title(args.proteinFile.split("/")[-1].split(".")[0]+"->"+args.target.split("/")[-1].split(".")[0])
    plt.xlabel("Normalized Score Values")
    plt.ylabel('Multiplicity')
    plt.savefig(args.outfiles+".pdf")



    outf=open(args.outfiles+".txt","w")
    json.dump(forHist,outf)
    
    outf.close()
   
    if not betterFile == None:
        while not q.empty():
            betterFile.write(q.get())
        betterFile.close()











