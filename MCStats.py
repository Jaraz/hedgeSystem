# -*- coding: utf-8 -*-
"""
Created on Sat Jan 10 21:49:14 2015

@author: jaycw_000
"""


class MCstats:
    def __init__(self):
        self.statVec = []
        self.deltaVec = []
    
    def addStats(self,val):
        self.statVec.append(val)
    
    def addDelta(self,val):
        self.deltaVec.append(val)
    
    def returnPrice(self):
        return sum(self.statVec) / len(self.statVec)
    
    def returnDelta(self):
        return sum(self.deltaVec) / len(self.deltaVec)
    
    def returnConvergence(self):
        n = len(self.statVec)
        answer = []
        
        while n > 10000:
            answer.append(sum(self.statVec[0:n]) / n) 
            n = n // 2

        answer.reverse()        
        
        return answer

    def returnStdError(self):
        n = len(self.statVec)
        answer = []
        
        while n > 100:
            answer.append(numpy.std(self.statVec[0:n])/n)
            n = n // 2
        
        answer.reverse()
        
        return answer