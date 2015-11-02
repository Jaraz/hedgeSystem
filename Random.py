# -*- coding: utf-8 -*-
"""
Created on Sat Jan 10 21:45:53 2015

@author: jaycw_000
"""
from __future__ import division
import numpy as np
from IPython import get_ipython
ipython = get_ipython()
import scipy.stats
import math

class randMC:
    def __init__(self, antiFlag, momentFlag) :
        self.antithetic = antiFlag
        self.momentMatching = momentFlag
        
    def genNormal(self, num):
        array = np.zeros(num)
        loopNum = num        
        
        #Check antithetic is even number
        if num % 2 == 0 and self.antithetic == True:
            loopNum = num // 2
        
        array[0:loopNum] = np.random.randn(loopNum)

        if num % 2 == 0 and self.antithetic == True:
            array[loopNum:num] = - array[0:loopNum]

        return array
        
    def genNormalMatrix(self, num, x, y):
        array = self.genNormal(num)        
        
        answer = np.zeros([x, y])
        
        for i in xrange(x):
            answer[i] = array[(i*y):((i+1)*y)]
        
        return answer        

class mcCDF:
    def __init__(self, cdfType):
        self.cdfType = cdfType
    
    def cdf(self,x):
        if self.cdfType == "approx":
            return normCDF(x)
        elif self.cdfType == "scipy":
            return scipy.stats.norm.cdf(x)
    

def normCDF(x):
    if x > 0:
        k = 1 / (1 + 0.2316419 * x)
        return 1 - 1 / (2 * math.pi ) ** 0.5 * math.exp( - x**2 / 2) * k * (
            0.319381530 + k * (
            -0.356563782 + k *(
            1.781477937 + k * (
            -1.821255978 + 1.330274429 * k))))
    if x < 0:
        return 1 - normCDF(-x)
    if x == 0:
        return 0.5

#vecCDF = np.vectorize(normCDF)
vecCDF = np.vectorize(scipy.stats.norm.cdf)