# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 20:18:52 2015

@author: jaycw_000
"""

from __future__ import division
import math
from scipy.stats import norm
from matplotlib.pyplot import *

class bsEngine:
    def __init__(self, antiFlag, momentFlag):
        self.rnd = randMC(antiFlag, momentFlag)
    
    def bsEvo(self, S, r, vol, optionType, num, stat):
        rndNumbers = self.rnd.genNormal(num)
        answer = 0.0    
        T = optionType.returnExpiry()
        discounting = math.exp(-r * T)
        term1 = r * T
        term2 = 0.5 * vol * vol * T
        term3 = vol * math.sqrt(T)
    
        for i in xrange(num):    
            term = math.exp(term1 - term2 + term3 * rndNumbers[i])
            S_final = S * term
            thisPayoff = optionType.payoff(S_final)
            answer += thisPayoff
            stat.addStats(discounting * thisPayoff)
            
            deltaAnswer = thisPayoff * (rndNumbers[i] / (S * term3))
            stat.addDelta(discounting * deltaAnswer)
              
        return discounting * answer / num


class eulerEngine:
    def __init__(self, antiFlag, momentFlag):
        self.rnd = randMC(antiFlag, momentFlag)

    def eulerEvo(self, S, r, vol, optionType, numSteps, numSims, stat):
        T = optionType.returnExpiry()
        delta_t = T / numSteps    
        discounting = math.exp(-r*T)
        answer = 0    
        rndNum = [[] for x in xrange(numSims)]
        volDeltaT = vol * math.sqrt(delta_t)
        rndNumbers = self.rnd.genNormal(numSteps * numSims)
        
        for i in xrange(numSims):
            for j in xrange(numSteps):
                rndNum[i].append(rndNumbers[i*numSteps + j])

        for i in xrange(numSims):
            s_next = 0
            s_last = S
        
            for j in xrange(numSteps):
                s_next = s_last + r * s_last * delta_t + s_last * volDeltaT * rndNum[i][j]
                s_last = s_next
            
            thisPayoff = optionType.payoff(s_next)
            stat.addStats(discounting * thisPayoff)
            answer += thisPayoff
    
        return discounting * answer / numSims