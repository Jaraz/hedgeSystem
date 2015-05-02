# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 20:18:52 2015

@author: jaycw_000
"""

from __future__ import division
import math
import numpy
import Random
from scipy.stats import norm
from scipy import integrate
from matplotlib.pyplot import *
from securities import *

class bsEngine:
    def __init__(self, antiFlag, momentFlag):
        self.rnd = Random.randMC(antiFlag, momentFlag)
    
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

class bachEngine:
    def __init__(self, antiFlag, momentFlag):
        self.rnd = Random.randMC(antiFlag, momentFlag)
    
    def bachEvo(self, S, r, vol, optionType, sims):
        rndNumbers = self.rnd.genNormal(sims)
        answer = numpy.zeros(sims)
        T = optionType.returnExpiry()
        
        answer = S * numpy.exp(r * T) + vol * math.sqrt(T) * rndNumbers
        
        return optionType.payoff(answer).mean()
        
    def bachNorm(self, S, r, vol, optionType, sims):
        rndNumbers = self.rnd.genNormal(sims)
        answer = numpy.zeros(sims)
        T = optionType.returnExpiry()
        
        answer = S * numpy.exp(r * T) + vol * numpy.sqrt((numpy.exp(2 * r * T) - 1)/(2*r)) * rndNumbers
        
        return optionType.payoff(answer).mean()
    
    def eulerStep(self, S, r, vol, optionType, steps, sims):
        rndNumbers = self.rnd.genNormalMatrix(steps * sims, sims, steps)
        rndNumbers = numpy.transpose(rndNumbers)
        answer = numpy.zeros(sims)
        answer[:] = S
        T = optionType.returnExpiry()
        dt = T / steps
        rdt = r * dt        
        rndVolDT = vol * numpy.sqrt(dt) * rndNumbers
                
        for i in xrange(steps):
            answer = answer + answer * rdt + rndVolDT[i]
        
        return optionType.payoff(answer).mean()
        

class bsIntegralEngine:
    def bsPrice(self, S, r, vol, optionType):
        T = optionType.returnExpiry()
        
        def integrand(x, S, expiry, vol, optionType, r):
            d = (numpy.log(x/S) + r * expiry + 0.5 * vol**2 * expiry) / (vol * numpy.sqrt(expiry))
            
            return optionType.payoff(x) / (vol * numpy.sqrt(expiry) * x) * norm.pdf(d)
        
        return(numpy.exp(-r * T) * integrate.quad(integrand, a = 0, b = 400, args = (S, T, vol, optionType, r), limit = 200)[0])

class eulerEngine:
    def __init__(self, antiFlag, momentFlag):
        self.rnd = Random.randMC(antiFlag, momentFlag)

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
        
        
## Main test Bach Model ##
#callOption = call(100, 1, 1)
callOption = call(100, 1, 1)
normEngine = bachEngine(True, False)
blackEngien= bsEngine(True, False)
port = portfolio([callOption])
bsInt = bsIntegralEngine()
resultStats = MCstats()

print "Analytic Call Norm Price      = ", bachPrice(1, 100, "call", 100, 0.1, 10)
print "MC Simulation Norm Call Price = ", normEngine.bachEvo(100, 0.1, 10, callOption, 1000000)
#print "MC 2                          = ", normEngine.bachNorm(100, 0.00000001, 10, callOption, 1000000)
print "MC Euleter Step Norm Call     = ", normEngine.eulerStep(100, 0.1, 10, callOption, 200, 1000000)
#print "MC Sim BS Price               = ", blackEngien.bsEvo(100, 0, .1, callOption, 100000, resultStats)
print "Analytic Call BS Price        = ", calcPrice(100, 0.1, .1, port)
#print "Integral Call BS Price        = ", bsInt.bsPrice(100, 0.01, .1, callOption)