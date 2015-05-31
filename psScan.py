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
import sys , os, ppr_phmm
'''
a wrapper for  PSscan. It also massages motif annotations and build s(6,1') for alignment. 
'''


def runScan(inp_):

    print"Searching For Motifs:"
    inp = inp_

    inputProt=open(inp,"r")

    prot=""


    #get seq info
    for line in inputProt:
        if not ">" in line:
            prot+=line.rstrip()
        
        
    #####Note to user:  hardcoded path change this to wherever ps_scan is ####################################################################################
    os.chdir("/s/parsons/b/others/approve/aPPRove/ps_scan")

    #run subprocess
    os.system("perl ps_scan.pl -d pprPrf.dat -l -1 -v -w pfsearch -C 250  -o pff "+inp+" > ../read.out")
    os.chdir("..")
    f = open("read.out","r")
    lines=[]
    for line in f:
        lines.append(line.rstrip())

    # get motif starting positions
    starts =[]

    for line in lines:
        if not ">" in line:
            vec=line.split("\t")
            if not prot[int(vec[1])-1]=="W":
                starts.append((int(vec[1])-1,int(vec[2])-1))
            else:
                starts.append((int(vec[1])-5,int(vec[2])-1))    
    starts.sort()
    
    
    #got to be valid motifs
    for i in range(len(starts)-1):
        if starts[i+1][0]-starts[i][1] > 30:
            starts.append((starts[i][1]+1,starts[i+1][0]-1))
            starts.sort()
             
           
    #print out start positions and get S6 and S1'
    print "Motif Locations:"    
    print starts
    ones=[]
    sixes=[]
    for startPos in starts:
        st = startPos
        
        ones.append(prot[st[0]])
        sixes.append(prot[st[0]+5])
    

    #The ending domains are also annotated but do not bind
    ones.pop(len(ones)-1)
    ones.pop(len(ones)-1)
    sixes.pop(len(sixes)-1)
    sixes.pop(len(sixes)-1)
    
    while((not ones[-1]=="D") and (not ones[-1]=="N") and (not ones[-1]=="S") ):
        ones.pop(len(ones)-1)
        sixes.pop(len(sixes)-1)

    motifs = []
    for i in range(len(sixes)-1):
        motifs.append(ones[i+1]+":"+sixes[i])
    
    types = [] 
    #assign types based off of size and orientation and were good to go   
    for i in range(len(sixes)-1):
        
        if (starts[i+1][0]-1-starts[i][0]) <= 33:
            types.append("S")
            
                    
            
        
        else:
            types.append("P")
            
                
    for i in range(len(types)):
       if types[i] == "S" and types[i-1] == "P" and types[i-2] == "P":
           types[i-1] = "L"                   
    
    print types
    
        
    for i in range(len(motifs)):
        motifs[i]= motifs[i]+types[i]
                
    print "Motif Binding Pairs:" 
       
    print motifs
    return motifs
    






'''
print "MEF1 rsp4" 
temp = "GTTGACTATGAAGAGAAGAATCAAAAGGATCGAACTACCTACTCATTATT"
ppr_phmm.run(temp,["I:C","S:N","D:S","D:P","N:S","N:S","S:S","N:S","T:L"])

print "MEF11 cox3"
temp ="TGAGGTTTTAGATCCTTGGGAAATCCCTTTTCTTAATACCCCTATTCTCC"
ppr_phmm.run(temp,["P:G","N:N","N:T","D:P","N:C","N:N","D:C","N:N","S:C","T:S","N:S","N:N","G:V"])
'''

if __name__ == "__main__":
    runScan(sys.argv[1])
    

#end of program    
#os.system("rm read.txt")
