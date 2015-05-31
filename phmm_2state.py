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
import sys,numpy
'''
This is the pair hidden markov model.
It is the algorithmic workhorse of aPPRove.


Thomas Harrison
3/19/2014
'''
class Node:
    """
    Container for dpMatrix used in the viterbi algorithm
    """
    def __init__(self,score=-1,state="m",ptr=None):
        # alignment score
        self.score=score
        # matrix symbol
        self.state=state
        # pointer to prev
        self.ptr=ptr
        ###################
        self.matchesCount=0
        ###################
        self.matchScore=0
        self.subAlignmentCount=0
        self.scorefromMatches=0
        
    def __repr__(self):
        return str(self.score)
        



def viterbi(s1,s2,emissions,emissionsGap,nu,tau,delta,epsilon,_z):
    '''
    s1 and s2 = sequence strings with s2<=s1 (for a two state PHMM)
    emissions = dictionary of binding probabilities ex key : AA with val .2
    delta = the probability from going from a match state to a gap state
    epsilon = the probability of a gap state to transfer back to a gap state
    gapCost = the probability of a gap in the alignment
    '''
    
    
    n = len(s1)
    m = len(s2)
    z=_z
    
    #initialize vM
    vM=[]
    for i in range (n+1):
        vM.append([])
        for j in range (m+1):
            vM[i].append([Node(float('-inf'),"m",None)])
            
            
    
    
    
    #vX
    vX=[]
    for i in range (n+1):
        vX.append([])
        for j in range (m+1):
            vX[i].append([Node(float('-inf'),"x",None)])
            
    
    #vD
    vD=[]
    for i in range (n+1):
        vD.append([])
        for j in range (m+1):
            vD[i].append([Node(float('-inf'),"x",None)])
            
    vDTwo=[]
    for i in range (n+1):
        vDTwo.append([])
        for j in range (m+1):
            vDTwo[i].append([Node(float('-inf'),"x",None)])            
    
    
    
    vD[0][0]=[Node(1.0,"x",None)]
    
    prev=vD[0][0][0]
    
    
    
    
    
    
    for i in range(1,n):
        
        vD[i][0]=[Node(prev.score+numpy.log(1-nu)+numpy.log(emissionsGap[s1[i-1]]),"d",prev)] 
        prev = vD[i][0][0]
            
        
    #vX[0][0][0].score=0.0
        
        
    
     
 
    for i in range(1,(n+1)) :
        for j in range(1,(m+1)):
        
            
         

 
 
        
            #max for transfering to a match state in vM
            prob = -1
            if(s1[i-1]+s2[j-1]) in emissions:
                prob = emissions[s1[i-1]+s2[j-1]]
            else:
                prob = emissions[s2[j-1]+s1[i-1]]
                
            possibleMoves=[]
            
            for node in vM[i-1][j-1]:
             
                ##################
                toAdd = numpy.log(1-delta-tau-0.0000001)+numpy.log(prob)
                ###################
                #Note for standardization
                #toAdd = numpy.log(1-delta-tau)+numpy.log(prob)
                         
                possibleMoves.append(Node((node.score)+toAdd,"m",node))
                possibleMoves[-1].matchesCount=node.matchesCount+1
                possibleMoves[-1].matchScore=(node.matchScore)+toAdd
                possibleMoves[-1].subAlignmentCount=node.subAlignmentCount+1
                
                
            if j==1:    
                for node in vD[i-1][j-1]:
                      
                    toAdd = numpy.log(nu)+numpy.log(prob)
                    possibleMoves.append(Node((node.score)+toAdd,"m",node))
                    possibleMoves[-1].matchesCount=node.matchesCount
                    possibleMoves[-1].matchScore=(node.matchScore)+toAdd-numpy.log(nu)
                    possibleMoves[-1].subAlignmentCount=node.subAlignmentCount                
                    
            for node in vX[i-1][j-1]:
                toAdd= numpy.log(1-epsilon)+numpy.log(prob)
                possibleMoves.append(Node((node.score)+toAdd,"m",node))
                possibleMoves[-1].matchesCount=node.matchesCount+1
                possibleMoves[-1].matchScore=(node.matchScore)+toAdd
                possibleMoves[-1].subAlignmentCount=node.subAlignmentCount+1                
               
     
                
            possibleMoves= sorted(possibleMoves, key=lambda node: node.score, reverse=True)
            vM[i][j]=[]
            k=0
            while k < len(possibleMoves) and k< z:
                vM[i][j].append(possibleMoves[k])
                k+=1
          
            #vDTwo
            possibleMoves=[]    
            
                  
            for node in vM[i-1][j]:
              
                toAdd = numpy.log(tau)+numpy.log(emissionsGap[s1[i-1]])
                    
                  
                possibleMoves.append(Node((node.score)+toAdd,"d",node))
                possibleMoves[-1].matchesCount=node.matchesCount
                possibleMoves[-1].matchScore=(node.matchScore)
                possibleMoves[-1].subAlignmentCount=node.subAlignmentCount
                    
            for node in vD[i-1][j]:
                ###############
                toAdd=numpy.log(1-nu)+numpy.log(emissionsGap[s1[i-1]])
                ########################
                #toAdd=numpy.log(1)+numpy.log(emissionsGap[s1[i-1]])
                possibleMoves.append(Node((node.score)+toAdd,"d",node))
                possibleMoves[-1].matchesCount=node.matchesCount                   
                possibleMoves[-1].matchScore=(node.matchScore)
                possibleMoves[-1].subAlignmentCount=node.subAlignmentCount                                     

            possibleMoves= sorted(possibleMoves, key=lambda node: node.score, reverse=True)
          
            vD[i][j]=[]
            k=0
            while k < len(possibleMoves) and k< z :
                vD[i][j].append(possibleMoves[k])
                k+=1
                
            possibleMoves=[]
            
            
                
            
            #vX      
            for node in vM[i-1][j]:
                toAdd = numpy.log(delta)+numpy.log(emissionsGap[s1[i-1]])
                possibleMoves.append(Node((node.score)+toAdd,"x",node))
                possibleMoves[-1].matchesCount=node.matchesCount 
                possibleMoves[-1].matchScore=(node.matchScore)+toAdd
                possibleMoves[-1].subAlignmentCount=node.subAlignmentCount+1                 
                    
            for node in vX[i-1][j]:
                toAdd=numpy.log(epsilon)+numpy.log(emissionsGap[s1[i-1]])
                possibleMoves.append(Node((node.score)+toAdd,"x",node))
                possibleMoves[-1].matchesCount=node.matchesCount
                possibleMoves[-1].matchScore=(node.matchScore)+toAdd
                possibleMoves[-1].subAlignmentCount=node.subAlignmentCount+1                                                       

            possibleMoves= sorted(possibleMoves, key=lambda node: node.score, reverse=True)
            vX[i][j]=[]
            k=0
            while k < len(possibleMoves) and k< z :
                vX[i][j].append(possibleMoves[k])
                k+=1
      
    res=[]                
    k=0
    
    
    #mToString(vM)
    #mToString(vD)
     
    for k in range(z):
        if k >= len(vM[-1][-1]):
            break
        #print mToString(vX)
        #print mToString(vM)
        ############
        vD[-1][-1][k].score+=numpy.log(nu)
        vM[-1][-1][k].score+=numpy.log(0.0000001)
        #############
        maxNode = vD[-1][-1][k]
        
        if vD[-1][-1][k].score < vM[-1][-1][k].score:
            maxNode = vM[-1][-1][k]
        
        #with division    
        res.append((vTrace(vM,vD,s1,s2,k),maxNode.matchScore/maxNode.subAlignmentCount,maxNode.score))
        #without division
        #res.append((vTrace(vM,vD,s1,s2,k),maxNode.matchScore,maxNode.score))
        
    return res
             
                  
             
            
            
                     
               
    
    
    
                 
     
                
     
def vTrace(vM,vD,s1,s2,k):

    #keeps track of the states through the optimal path
    stateStr = ""
    #first grab the max of the nth row and mth column of both matricies
    cur = None
    if vM[-1][-1][k].score > vD[-1][-1][k].score:
        cur = vM[-1][-1][k]
    else:
        cur = vD[-1][-1][k]

    
    #Trace the pointers back until it reaches the terminal point and append to list
    while(not cur.ptr == None):
        stateStr+=cur.state
        cur = cur.ptr
    #flip the state string
    stateStr = stateStr[::-1]

    #given the state string make a string representative of the alignment
   

    al=""
    if len(s2[0])== 1:
        alOne=""
        alTwo=""
        al=""
        countOne=0
        countTwo=0
        for char in stateStr:
            if char == "m":
                alOne+=s1[countOne]
                alTwo+=s2[countTwo]
                countOne+=1
                countTwo+=1
            else:
                alOne+=s1[countOne]
                countOne+=1
                alTwo+="_"
        al= alOne+"\n"+alTwo
        
    else:
        onePrime="1'"
        six="6 "
        N = "T "
        al=""
        countOne=0
        countTwo=0
        for char in stateStr:
            if char == "m":
                tokens = s2[countTwo].split(":")
                onePrime+=tokens[0]
                six+= tokens[1][0]
                N+=s1[countOne]
                countOne+=1
                countTwo+=1
            else:
                if len(stateStr)<= 50:
                    N+=s1[countOne]                
                countOne+=1
                if len(stateStr)<=50:
                    onePrime+="_"
                    six+="_"
        al= onePrime+"\n"+six+"\n"+N
     
     
    #print stateStr
    #return the alignment and the state string
    return (al , stateStr)
     
           

def mToString(m):
    st = ""
    
    for item in m:
        for test in item:
             st+= str(test[0])+str(test[0].state)+","
        st += "\n"
    print st


#main for debugging
if __name__ == "__main__":

   
    

    
    emissions = {"AA": .2, "AG": .2/4, "AC": .2/4, "AT": .2/4, "GG": .2, "GC": .2/4 , "GT":.2/4, "CC": .2, "CT" :.2/4, "TT":.2}
    gapCost = .2/4
    delta = .5
    epsilon = .5

    s1= "AGCTCAAGTCATCA"
    s2= "CAAG"
    print "Test: "+s1+" and "+s2+"."
    res = viterbi(s1,s2,emissions,delta,epsilon,gapCost,3)
    
    print "Best Al"
    print res[0][0][0]
    print "with Score: "+str(res[0][1])
    
    print "Second Best Al"
    print res[1][0][0]
    print "with Score: "+str(res[1][1])
    
    print "Third Best Al"
    print res[2][0][0]
    print "with Score: "+str(res[2][1])
    
    
    s1= "AGCTCAAGTCATCA"
    s2= "AGC"
    print "Test: "+s1+" and "+s2+"."
    res = viterbi(s1,s2,emissions,delta,epsilon,gapCost,20)
    
    print "Best Al"
    print res[0][0][0]
    print "with Score: "+str(res[0][1])
    
    print "Second Best Al"
    print res[1][0][0]
    print "with Score: "+str(res[1][1])
    
    print "Third Best Al"
    print res[2][0][0]
    print "with Score: "+str(res[2][1])
    
    
    
    '''
    print "Alignment: \n"+al
    print "States: \n"+s
    print"__________________________\n"

    s1= "AGCTGA"
    s2= "TGA"
    print "Test: "+s1+" and "+s2+"."
    al,s=viterbi(s1,s2,emissions,delta,epsilon,gapCost)
    print "Alignment: \n"+al
    print "States: \n"+s
    print"__________________________\n"
    
    s1= "AGCTGA"
    s2= "CTG"
    print "Test: "+s1+" and "+s2+"."
    al,s=viterbi(s1,s2,emissions,delta,epsilon,gapCost)
    print "Alignment: \n"+al
    print "States: \n"+s
    print"__________________________\n"

    s1= "AGCTGA"
    s2= "ACG"
    print "Test: "+s1+" and "+s2+"."
    al,s=viterbi(s1,s2,emissions,delta,epsilon,gapCost)
    print "Alignment: \n"+al
    print "States: \n"+s
    print"__________________________\n"
    '''
     
                                  
            
     
    
    
    
    
    
    
    
    
