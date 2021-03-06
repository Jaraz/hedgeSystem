# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 09:47:36 2015

@author: jaycw_000
"""

import numpy

#Option Classes
class option:
    desc = ""    
    
    def __init__(self, K, T, x, european = True):
        self.expiry = T
        self.strike = K
        self.units = x
        self.european = european
    
    def updateExpiry(self,x):
        self.expiry = x
        
    def returnExpiry(self):
        return self.expiry
        
    def returnStrike(self):
        return self.strike
        
    def returnDesc(self):
        return self.desc
    
    def returnUnits(self):
        return self.units
        
    def returnEuropean(self):
        return self.european

class call(option): 
    desc = "call"
    
    def payoff(self, S):
        return self.units * numpy.maximum(S - self.strike, 0)

class put(option):
    desc = "put"
    
    def payoff(self, S):
        return self.units * numpy.maximum(self.strike - S, 0)

class callDigital(option):
    desc = "digCall"
    
    def payoff(self, S):
        answer = numpy.choose(S>self.strike, [0,1])
        return self.units * answer
        #if S > self.strike:
        #    return self.units * 1
        #else:
        #    return self.units * 0

class putDigital(option):
    desc = "digPut"
    
    def payoff(self, S):
        answer = numpy.choose(S<self.strike, [0,1])
        return self.units * answer
        #if S < self.strike:
        #    return self.units * 1
        #else:
        #    return self.units * 0

class forward(option):
    desc = "fwd"
    
    def payoff(self, S):
        return self.units * (S - self.strike)

class stock(option):
    desc = "stock"
    
    def payoff(self,S):
        return self.units * (S - 0)
        
#Portfolio
class portfolio:
    def __init__(self, optList):
        self.opt = []
        self.opt = optList

    def returnSec(self):
        return self.opt
        
    def addSec(self, addList):
        self.opt += addList
    
    def removeLast(self):
        self.opt = self.opt[0:len(self.opt)-1]
        
    def updateExpiry(self,x):
        for i in xrange(len(self.opt)):
            self.opt[i].updateExpiry(x)

    def returnExpiry(self):
        return self.opt[0].returnExpiry()

    def returnDesc(self):
        for i in xrange(len(self.opt)):
            print "Opt ", i, "desc: ", self.opt[i].returnDesc(), "exp: ", self.opt[i].returnExpiry(), "k: ", self.opt[i].returnStrike(), "unit: ", self.opt[i].returnUnits()