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
import phmm_2state as phmm
import random

'''
ppr-phmm messages the input sent in by approve for the phmm.
The parameterization of this model is based off barkan et al.

'''


emissions = {}
emissionsGap={}
nu=-1
tau=-1
delta = -1
epsilon = -1
z = 1
C = 0.0



def vTrain(d,gapCost):
    # in case viterbi training is needed
    global emissions
    global delta
    global epsilon

    
    totalMs=0.0
    fromMatch = 0.0
    fromMatchToGap=0.0
        
    fromGap=0.0
    fromGapToGap=0.0
    #zero out emissions dict
    emissions = dict.fromkeys(emissions, 0.0)
   
    
    for pair in d:
        
    
        
        

        
        s = pair[0]
        #get new transition probabilities
        for j in range(len(s)-1):
            if "m" in s[j]:
                fromMatch+=1
                if "x" in s[j+1]:
                    fromMatchToGap +=1
            if "x" in s[j]:
                fromGap+=1
                if "x" in s[j+1]:
                    fromGapToGap +=1
           

                   
                   
        #get new emission probabilities
        
        #get all 3 (1', 6, and target)   
        alPieces = pair[1]
        
        
        
        for i in range(2,len(alPieces[0])):
         
            if not "_" in alPieces[0][i]:
       
        
                #Get key of binding and add 1    
                emissions[alPieces[0][i]+":"+alPieces[1][i]+alPieces[2][i]] +=1
                totalMs += 1
                
                
    delta = fromMatchToGap/fromMatch       
    epsilon = fromGapToGap/fromGap            
    #divide all emission probabilities by the total number of emissions.          
    for item in emissions:
        emissions[item] /= totalMs
    
    
        
        
    
       
            
        
    
def setParams(_nu,_tau,_delta,_epsilon,_z):
    # parametrizes the model of phmm
    global nu
    global tau
    global delta
    global epsilon
    global emissions
    global emissionsGap
    global z
    global C
    nu=_nu
    tau=_tau
    delta=_delta
    epsilon=_epsilon
    z=_z
    aaString = "ACDEFGHIKLMNPQRSTVWY"
    bpString = "AUGC"    
    #Create emmision keys for dict
    for i in range(len(aaString)):
        for j in range(len(aaString)):
            for k in range(len(bpString)):
                key = aaString[i]+":"+aaString[j]+"P"+bpString[k]
                emissions[key]=0.0/1600
                key = aaString[i]+":"+aaString[j]+"L"+bpString[k]
                emissions[key]=0.0/1600
                key = aaString[i]+":"+aaString[j]+"S"+bpString[k]
                emissions[key]=0.0/1600
                
                   
                
    #PMotifs                         
    emissions["D:TPA"]+=2 
    emissions["D:TPG"]+=23 
    emissions["D:TPU"]+=1
    emissions["N:TPA"]+=23
    emissions["N:TPG"]+=2
    emissions["D:NPA"]+=7
    emissions["D:NPC"]+=22 
    emissions["D:NPG"]+=7 
    emissions["D:NPU"]+=68
    emissions["N:NPA"]+=6 
    emissions["N:NPC"]+=46 
    emissions["N:NPG"]+=3 
    emissions["N:NPU"]+=31
    emissions["N:SPA"]+=10 
    emissions["N:SPC"]+=1
    emissions["N:SPU"]+=1
    emissions["S:NPA"]+=1 
    emissions["S:NPC"]+=12 
    emissions["S:NPU"]+=3
    emissions["T:TPG"]+=3
    emissions["S:SPA"]+=5 
    emissions["S:SPC"]+=1
    emissions["R:SPG"]+=2
    emissions["C:SPA"]+=2
    emissions["D:RPU"]+=3
    emissions["S:TPA"]+=3
    emissions["S:TPU"]+=1
    emissions["R:NPC"]+=2
    emissions["D:APG"]+=2
    emissions["D:APU"]+=1
    emissions["D:GPC"]+=1 
    emissions["D:GPG"]+=2
    emissions["M:HPU"]+=2
    emissions["S:HPU"]+=2
    emissions["C:NPC"]+=5 
    emissions["C:NPU"]+=3
    emissions["G:FPA"]+=1
    emissions["H:GPA"]+=1
    emissions["S:GPA"]+=1
    emissions["G:SPA"]+=1
    emissions["L:SPA"]+=1
    emissions["P:TPA"]+=1
    emissions["V:HPA"]+=2 
    emissions["V:HPU"]+=1
    emissions["N:GPC"]+=1
    emissions["T:NPA"]+=1
    emissions["T:NPC"]+=8 
    emissions["T:NPU"]+=8
    emissions["D:NPA"]+=1
    emissions["D:NPC"]+=2 
    emissions["D:NPU"]+=5
    emissions["N:APA"]+=2 
    emissions["N:APC"]+=2
    emissions["D:MPC"]+=1 
    emissions["D:MPG"]+=1
    emissions["E:NPG"]+=1 
    emissions["E:NPU"]+=1
    emissions["R:TPA"]+=1 
    emissions["R:TPG"]+=1
    emissions["D:SPA"]+=4
    emissions["D:SPC"]+=3 
    emissions["D:SPG"]+=3 
    emissions["D:SPU"]+=1
    emissions["S:CPC"]+=1 
    emissions["S:CPU"]+=2
    emissions["D:IPA"]+=1 
    emissions["D:IPG"]+=1 
    emissions["D:IPU"]+=1
        
    #SMotifs
    emissions["D:TSA"]+=9
    emissions["D:TSC"]+=2 
    emissions["D:TSG"]+=27 
    emissions["D:TSU"]+=4
    emissions["D:SSA"]+=1
    emissions["D:SSC"]+=1
    emissions["D:SSG"]+=13
    emissions["D:SSU"]+=1
    emissions["N:SSA"]+=20
    emissions["N:SSC"]+=2
    emissions["N:SSG"]+=1
    emissions["N:SSU"]+=5
    emissions["N:TSA"]+=18
    emissions["N:TSC"]+=1
    emissions["N:TSG"]+=2
    emissions["N:TSU"]+=3
    emissions["D:NSA"]+=5
    emissions["D:NSC"]+=13
    emissions["D:NSG"]+=5
    emissions["D:NSU"]+=44
    emissions["N:CSA"]+=3
    emissions["T:NSC"]+=12
    emissions["T:NSG"]+=3
    emissions["T:NSU"]+=4
    emissions["N:NSA"]+=5
    emissions["N:NSC"]+=8
    emissions["N:NSG"]+=9
    emissions["N:NSU"]+=5
    emissions["S:GSA"]+=2
    emissions["R:ASC"]+=2
    emissions["H:TSA"]+=1
    emissions["H:TSG"]+=2
    emissions["D:CSG"]+=1
    emissions["D:GSG"]+=1
    emissions["D:GSG"]+=1
    emissions["C:PSG"]+=1
    emissions["H:ASA"]+=1
    emissions["E:SSA"]+=1
    emissions["E:TSA"]+=1
    emissions["K:TSA"]+=1
    emissions["K:SSA"]+=2
    emissions["K:SSG"]+=1
    emissions["L:TSA"]+=2
    emissions["L:TSU"]+=1
    emissions["D:LSC"]+=1
    emissions["E:NSC"]+=1
    emissions["T:ASU"]+=1
    emissions["N:FSU"]+=1
    emissions["N:ISU"]+=1
    emissions["N:PSU"]+=1
    emissions["S:SSU"]+=1
    emissions["Q:TSU"]+=1
    emissions["S:VSU"]+=1
    emissions["D:KSC"]+=2
    emissions["D:KSG"]+=1
    emissions["H:NSC"]+=2
    emissions["H:NSU"]+=1
    emissions["D:ASG"]+=1
    emissions["D:ASU"]+=1
    emissions["S:TSA"]+=1
    emissions["S:TSC"]+=3
    emissions["S:TSG"]+=1
    emissions["S:TSU"]+=1
    emissions["N:ASA"]+=1
    emissions["N:ASU"]+=1
    emissions["T:TSA"]+=1
    emissions["T:TSC"]+=1
    emissions["T:TSG"]+=1
    emissions["S:NSA"]+=2
    emissions["S:NSC"]+=4
    emissions["S:NSG"]+=1
    emissions["S:NSU"]+=5
        
        
    #get probs with a psuedocount of one
    for token in emissions:
        emissions[token]+=1
        C+=emissions[token]
            
    for token in emissions:
        emissions[token] /= C
        
    total = 0.0
    
    #gap emmisions should all be the same

    emissionsGap["A"] = 0.25
    emissionsGap["C"] = 0.25
    emissionsGap["G"] = 0.25
    emissionsGap["U"] = 0.25
    
 
def filterOut(filterList):
    # in case a user want LOO cross validation
    global C
    global emssions
    C=0
    aaString = "ACDEFGHIKLMNPQRSTVWY"
    bpString = "AUGC"    
    #Create keys for dict
    for i in range(len(aaString)):
        for j in range(len(aaString)):
            for k in range(len(bpString)):
                key = aaString[i]+":"+aaString[j]+"P"+bpString[k]
                emissions[key]=0.0/1600
                key = aaString[i]+":"+aaString[j]+"L"+bpString[k]
                emissions[key]=0.0/1600
                key = aaString[i]+":"+aaString[j]+"S"+bpString[k]
                emissions[key]=0.0/1600
                
                   
                
    #PMotifs                         
    emissions["D:TPA"]+=2 
    emissions["D:TPG"]+=23 
    emissions["D:TPU"]+=1
    emissions["N:TPA"]+=23
    emissions["N:TPG"]+=2
    emissions["D:NPA"]+=7
    emissions["D:NPC"]+=22 
    emissions["D:NPG"]+=7 
    emissions["D:NPU"]+=68
    emissions["N:NPA"]+=6 
    emissions["N:NPC"]+=46 
    emissions["N:NPG"]+=3 
    emissions["N:NPU"]+=31
    emissions["N:SPA"]+=10 
    emissions["N:SPC"]+=1
    emissions["N:SPU"]+=1
    emissions["S:NPA"]+=1 
    emissions["S:NPC"]+=12 
    emissions["S:NPU"]+=3
    emissions["T:TPG"]+=3
    emissions["S:SPA"]+=5 
    emissions["S:SPC"]+=1
    emissions["R:SPG"]+=2
    emissions["C:SPA"]+=2
    emissions["D:RPU"]+=3
    emissions["S:TPA"]+=3
    emissions["S:TPU"]+=1
    emissions["R:NPC"]+=2
    emissions["D:APG"]+=2
    emissions["D:APU"]+=1
    emissions["D:GPC"]+=1 
    emissions["D:GPG"]+=2
    emissions["M:HPU"]+=2
    emissions["S:HPU"]+=2
    emissions["C:NPC"]+=5 
    emissions["C:NPU"]+=3
    emissions["G:FPA"]+=1
    emissions["H:GPA"]+=1
    emissions["S:GPA"]+=1
    emissions["G:SPA"]+=1
    emissions["L:SPA"]+=1
    emissions["P:TPA"]+=1
    emissions["V:HPA"]+=2 
    emissions["V:HPU"]+=1
    emissions["N:GPC"]+=1
    emissions["T:NPA"]+=1
    emissions["T:NPC"]+=8 
    emissions["T:NPU"]+=8
    emissions["D:NPA"]+=1
    emissions["D:NPC"]+=2 
    emissions["D:NPU"]+=5
    emissions["N:APA"]+=2 
    emissions["N:APC"]+=2
    emissions["D:MPC"]+=1 
    emissions["D:MPG"]+=1
    emissions["E:NPG"]+=1 
    emissions["E:NPU"]+=1
    emissions["R:TPA"]+=1 
    emissions["R:TPG"]+=1
    emissions["D:SPA"]+=4
    emissions["D:SPC"]+=3 
    emissions["D:SPG"]+=3 
    emissions["D:SPU"]+=1
    emissions["S:CPC"]+=1 
    emissions["S:CPU"]+=2
    emissions["D:IPA"]+=1 
    emissions["D:IPG"]+=1 
    emissions["D:IPU"]+=1
        
    #SMotifs
    emissions["D:TSA"]+=9
    emissions["D:TSC"]+=2 
    emissions["D:TSG"]+=27 
    emissions["D:TSU"]+=4
    emissions["D:SSA"]+=1
    emissions["D:SSC"]+=1
    emissions["D:SSG"]+=13
    emissions["D:SSU"]+=1
    emissions["N:SSA"]+=20
    emissions["N:SSC"]+=2
    emissions["N:SSG"]+=1
    emissions["N:SSU"]+=5
    emissions["N:TSA"]+=18
    emissions["N:TSC"]+=1
    emissions["N:TSG"]+=2
    emissions["N:TSU"]+=3
    emissions["D:NSA"]+=5
    emissions["D:NSC"]+=13
    emissions["D:NSG"]+=5
    emissions["D:NSU"]+=44
    emissions["N:CSA"]+=3
    emissions["T:NSC"]+=12
    emissions["T:NSG"]+=3
    emissions["T:NSU"]+=4
    emissions["N:NSA"]+=5
    emissions["N:NSC"]+=8
    emissions["N:NSG"]+=9
    emissions["N:NSU"]+=5
    emissions["S:GSA"]+=2
    emissions["R:ASC"]+=2
    emissions["H:TSA"]+=1
    emissions["H:TSG"]+=2
    emissions["D:CSG"]+=1
    emissions["D:GSG"]+=1
    emissions["D:GSG"]+=1
    emissions["C:PSG"]+=1
    emissions["H:ASA"]+=1
    emissions["E:SSA"]+=1
    emissions["E:TSA"]+=1
    emissions["K:TSA"]+=1
    emissions["K:SSA"]+=2
    emissions["K:SSG"]+=1
    emissions["L:TSA"]+=2
    emissions["L:TSU"]+=1
    emissions["D:LSC"]+=1
    emissions["E:NSC"]+=1
    emissions["T:ASU"]+=1
    emissions["N:FSU"]+=1
    emissions["N:ISU"]+=1
    emissions["N:PSU"]+=1
    emissions["S:SSU"]+=1
    emissions["Q:TSU"]+=1
    emissions["S:VSU"]+=1
    emissions["D:KSC"]+=2
    emissions["D:KSG"]+=1
    emissions["H:NSC"]+=2
    emissions["H:NSU"]+=1
    emissions["D:ASG"]+=1
    emissions["D:ASU"]+=1
    emissions["S:TSA"]+=1
    emissions["S:TSC"]+=3
    emissions["S:TSG"]+=1
    emissions["S:TSU"]+=1
    emissions["N:ASA"]+=1
    emissions["N:ASU"]+=1
    emissions["T:TSA"]+=1
    emissions["T:TSC"]+=1
    emissions["T:TSG"]+=1
    emissions["S:NSA"]+=2
    emissions["S:NSC"]+=4
    emissions["S:NSG"]+=1
    emissions["S:NSU"]+=5
        
    for key in filterList:
        if emissions[key]>=2:
            emissions[key]-=1
    #print emissions  
        
    
    for token in emissions:
        emissions[token]+=1
        C+=emissions[token]
            
    for token in emissions:
        emissions[token] /= C
        
    total = 0.0
    
    emissionsGap["A"] = 0.25
    emissionsGap["C"] = 0.25
    emissionsGap["G"] = 0.25
    emissionsGap["U"] = 0.25 
    


  
# run phmm
def run(temp,s2):
    global nu
    global tau
    global delta
    global epsilon
    global emissions
    global emissionsGap
    global z
    s1=""
    # convert T's to U's this is RNA after all...
    for char in temp:
        if "T" in char:
            s1+="U"
        else:
            s1+=char
    #T(D1,D1) ought to be dependent on transcript size.
    if not (len(s1)-len(s2))==0:
        nu = 1.0/((len(s1)-len(s2))/2.0)
    else:
        #hmmmm... no gaps.... well just give it a very small psuedocount then
        nu = 0.000000001
    #let it run    
    res=phmm.viterbi(s1,s2,emissions,emissionsGap,nu,tau,delta,epsilon,z)
    return res
    
    
   
    
        
        
    

