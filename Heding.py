# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 19:19:16 2015

@author: jaycw_000
"""

from __future__ import division
import math
from scipy.stats import norm
from matplotlib.pyplot import *
import numpy as np
import copy
from securities import *
from optFuncs import *
from Random import *
from IPython.parallel import Client

import pickle
import dill

from types import FunctionType
from IPython.utils.pickleutil import can_map
from IPython.kernel.zmq import serialize

can_map.pop(FunctionType, None)
serialize.pickle = pickle


def runHedgeSimul(S, T, port, r, vol, path, hedgeType, debugFlag):
    pnl = 0
    s_old = S
    s_new = S
    
    hedgeT = 0
    hedgeTm1 = 0

    numSteps = len(path)
    
    for i in xrange(0, numSteps):
        expiry = T * (1 - i/numSteps)
        port.updateExpiry(expiry)
        hedge = hedgeType.returnHedge(s_new, r, vol, port)
        
        #MTM of gamma hedge
        if hedge[1] <> 0:                
            atmHedge = call(s_new, expiry, hedge[1])
            atmPort = portfolio([atmHedge])
            hedgeTm1 = calcPrice(s_new, r, vol, atmPort)
        
        s_old = s_new
        s_new = path[i]
    
        if hedge[1] <> 0:
            if i < numSteps-1:        
                atmPort.updateExpiry(T * (1 - (i+1)/numSteps))
                hedgeT = calcPrice(s_new, r, vol, atmPort)            
            else:
                hedgeT = calcPayoff(s_new, atmPort)
         
        if debugFlag == True: 
            print "Step: ", i
            print "s t-1: ", s_old
            print "s t: ", s_new
            print "hedges: ", hedge
            print "hedgeTm1: ", hedgeTm1
            print "hedge: ", hedgeT
            print "stock pnl: ", hedge[0] * (s_new - s_old)
            print "gamma pnl: ", hedgeT - hedgeTm1
            print " "
                    
        #PnL from delta hedging
        pnl += hedge[0] * (s_new - s_old)
        
        #PnL from gamma hedging
        if hedge[1] <>0:
            pnl += hedgeT - hedgeTm1

    optionPayoff = calcPayoff(s_new, port)
    pnl = pnl + optionPayoff
    
    return pnl    

class gbmPathEngine:
    def __init__(self):
        self.rnd = randMC(True,False)

    def eulerPath(self, T, S, mu, r, vol, numSteps, numPaths):
        delta_t = T / numSteps
        
        answer = np.zeros([numPaths, numSteps])
        randNumbers = np.zeros([numPaths, numSteps])
        
        randNumbers = self.rnd.genNormalMatrix(numPaths * numSteps, numPaths, numSteps)
        
        volDeltaT = vol * math.sqrt(delta_t)        
        
        for j in xrange(numPaths):
            volDeltaRand = randNumbers[j] * volDeltaT
            
            term = 1.0 + mu * delta_t + volDeltaRand
    
            s_next = 0
            s_last = S
    
            for i in xrange(0,numSteps):
                s_next = s_last * term[i]
                s_last = s_next
                answer[j, i] = s_next
    
        return answer

#returns an array for gamma hedging
class hedger:
    def __init__(self, gammaFlag):
        self.gammaHedge = gammaFlag
    
    def returnHedge(self, S, r, vol, portfolio):
       
        if(S > K):
            return np.array([-1,0]);
        else:
            return np.array([0,0]);
   
class bsHedger(hedger):
    def returnHedge(self, S, r, vol, port):
        bsGamma = 0

        #newPort = copy.deepcopy(port)
        newPort = portfolio(port.returnSec())

        if self.gammaHedge == True:
            portGamma = calcGamma(S, r, vol, port)
            atmGamma = optionGamma(port.returnExpiry(), S, "call", S, r, vol)
            bsGamma = -portGamma / atmGamma        
            newSec = call(S, port.returnExpiry(), bsGamma)
            
            #newPort.addSec([newSec])
            newPort.opt = port.opt + [newSec]
        
        bsDelta = -calcDelta(S, r, vol, newPort)
        
        return np.array([bsDelta,bsGamma])

class hedgeSimul:
    def __init__(self, S, mu, r, vol, portfolio, debugFlag):
        self.S = S
        self.port = portfolio
        self.secList = portfolio.returnSec()
        self.T = self.secList[0].returnExpiry()
        self.mu = mu
        self.r = r
        self.vol = vol
        self.debugFlag = debugFlag
        self.hedgeType = hedger
        
        
    def setMu(self, x):
        self.mu = x
    
    def setVol(self,x):
        self.vol = x
    
    def runSim(self, numSteps, numSims, pathEngine, hedger):
        stats = []
        self.hedgeType = hedger
        path = pathEngine.eulerPath(self.T, self.S, self.mu, self.r, self.vol, numSteps, numSims)

        clients = Client()

        lview = clients.load_balanced_view()
        dview = clients[:]        
        
        with dview.sync_imports():
            import numpy as np
            from Random import randMC
            from securities import portfolio
            from securities import call
            from optFuncs import calcGamma
            from optFuncs import optionGamma
            from optFuncs import calcDelta
            from optFuncs import optionDelta
            from optFuncs import calcPrice
            from optFuncs import optionVal
        
        
        dview.push(dict(runHedgeSimul = runHedgeSimul))        
        
        myLambdaFunc = lambda path: runHedgeSimul(S, T, port, r, vol, path, hedgeType, debugFlag)

        S = self.S
        T = self.T
        r = self.r
        vol = self.vol
        port = self.port
        hedgeType = self.hedgeType
        debugFlag = self.debugFlag
        
        myDict = dict(S = S, T = T, r = r, vol = vol, port = port, hedgeType = hedgeType, debugFlag = debugFlag)
        dview.push(myDict)

        lview.block = True

        pickle.dumps(myLambdaFunc)
        

        stats = lview.map(myLambdaFunc, path)        
#runHedgeSimul(S, T, port, r, vol, path, hedgeType, debugFlag)        
        #for j in xrange(numSims):   
        #    stats.append(runHedgeSimul(self.S, self.T, self.port, self.r, self.vol, path[j], self.hedgeType, self.debugFlag))
            
        return stats
        
        
    def runPath(self, path):
        pnl = 0
        s_old = self.S
        s_new = self.S
        
        hedgeT = 0
        hedgeTm1 = 0
        
        numSteps = len(path)
        
        for i in xrange(0, numSteps):
            expiry = self.T * (1 - i/numSteps)
            self.port.updateExpiry(expiry)
            hedge = self.hedgeType.returnHedge(s_new, self.r, self.vol, self.port)
            
            #MTM of gamma hedge
            if hedge[1] <> 0:                
                atmHedge = call(s_new, expiry, hedge[1])
                atmPort = portfolio([atmHedge])
                hedgeTm1 = calcPrice(s_new, self.r, self.vol, atmPort)
            
            s_old = s_new
            s_new = path[i]
        
            if hedge[1] <> 0:
                if i < numSteps-1:        
                    atmPort.updateExpiry(self.T * (1 - (i+1)/numSteps))
                    hedgeT = calcPrice(s_new, self.r, self.vol, atmPort)            
                else:
                    hedgeT = calcPayoff(s_new, atmPort)
             
            if self.debugFlag == True: 
                print "Step: ", i
                print "s t-1: ", s_old
                print "s t: ", s_new
                print "hedges: ", hedge
                print "hedgeTm1: ", hedgeTm1
                print "hedge: ", hedgeT
                print "stock pnl: ", hedge[0] * (s_new - s_old)
                print "gamma pnl: ", hedgeT - hedgeTm1
                print " "
                        
            #PnL from delta hedging
            pnl += hedge[0] * (s_new - s_old)
            
            #PnL from gamma hedging
            if hedge[1] <>0:
                pnl += hedgeT - hedgeTm1

        optionPayoff = calcPayoff(s_new, self.port)
        pnl = pnl + optionPayoff
        
        return pnl


def conf_int_native(x, ci=0.95):
    ci2 = (1-ci)*.5
    low_idx = int(ci2*x.size)
    high_idx = int((1-ci2)*x.size)
    x.sort()
    return x.mean(), x[low_idx], x[high_idx]


#main
#clients = Client()
#
#lview = clients.load_balanced_view()
#dview = clients[:]
#
#with dview.sync_imports():
#    import numpy
#    from Random import randMC
#    from Hedging import * 
#
#
#
MCengine = gbmPathEngine()
hdgrDelta = bsHedger(False)
hdgrGamma = bsHedger(True)
callOption = call(110, 1, 100)
callATM = call(100, 1, 1)
numSims = 100

secList = [callOption]
port1 = portfolio(secList)

sim = hedgeSimul(100, 0.0, 0.0, 0.1, port1, False)
    
test = np.array(sim.runSim(40, numSims, MCengine, hdgrGamma))
print "PnL = ", sum(test)/numSims
print "sd = ", np.std(test)
print "conf = ", conf_int_native(test)

##hist(test)