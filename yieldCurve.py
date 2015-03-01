# -*- coding: utf-8 -*-
"""
Created on Sun Feb 15 15:35:26 2015

@author: jaycw_000
"""

from fredapi import Fred
import numpy
import scipy
import matplotlib
from matplotlib import pyplot
from scipy import interpolate
from scipy import optimize
fred = Fred(api_key='ba6be791f155502772efcac065904210')

class yieldCurve:
    def __init__(self, inputDate):
        self.curve = curveBuild(inputDate)
        temp = -numpy.log(self.curve[1])/self.curve[0]
        temp[0] = 0
        self.zeroCurve = interpolate.PchipInterpolator(self.curve[0], temp)
        
        
    def discFact(self, t):
        if t < 0 :
            return 1
        return numpy.exp(-self.zeroCurve(t)*t)

#pull swap data for a date
def pullSwapData(inputDate):
    swapArray = numpy.zeros(9)
    swapArray[0] = 0.2581
    swapArray[1] = fred.get_series('DSWP1', observation_start=inputDate, observation_end=inputDate)[0]
    swapArray[2] = fred.get_series('DSWP2', observation_start=inputDate, observation_end=inputDate)[0]
    swapArray[3] = fred.get_series('DSWP3', observation_start=inputDate, observation_end=inputDate)[0]
    swapArray[4] = fred.get_series('DSWP4', observation_start=inputDate, observation_end=inputDate)[0]
    swapArray[5] = fred.get_series('DSWP5', observation_start=inputDate, observation_end=inputDate)[0]
    swapArray[6] = fred.get_series('DSWP7', observation_start=inputDate, observation_end=inputDate)[0]
    swapArray[7] = fred.get_series('DSWP10', observation_start=inputDate, observation_end=inputDate)[0]
    swapArray[8] = fred.get_series('DSWP30', observation_start=inputDate, observation_end=inputDate)[0]
    
    return swapArray/100
    
def curveBuild(inputDate):
    data = pullSwapData(inputDate)
    curve = prepCurve()
    
    #use First deposit
    curve[1,1] = 1 / (1 + 0.25 * data[0])
    
    for j in xrange(3):
        for i in xrange(2,10):
            curve[1,i] = scipy.optimize.brentq(swapOptim, a = 0.01, b = 0.9999, args = (i, 0, curve[0,i], data[i-1], curve))
    
    return curve
    
#assume 3m libor and periods are always whole numbers
def swapPricer(start, tenor, strike, curve):
    fix_per = int(tenor * 2)
    flt_per = int(tenor * 4)
    
    flt = 0
    fix = 0
    
    for i in xrange(fix_per):    
        fix = fix + 0.5 * curveInterp(start + i * 0.5, curve) * strike
    
    for j in xrange(flt_per):
        flt = flt + 0.25 * curveInterp(start + j * 0.25, curve) * fwdRate(start + j * 0.25, curve)
    
    return fix - flt

def swapOptim(x, n, start, tenor, strike, curve):
    curve[1,n] = x
    
    return swapPricer(start, tenor, strike, curve)

def swapDV01(start, tenor, curve):
    fix_per = int(tenor * 2)
    
    dv01 = 0
    
    for i in xrange(fix_per):    
        dv01 = dv01 + 0.5 * curveInterp(start + i * 0.5, curve)
    
    return dv01

def prepCurve():
    x = numpy.zeros((2,10))
    x.fill(0.5)
    x[0,0] = 0
    x[0,1] = 0.25
    x[0,2] = 1
    x[0,3] = 2
    x[0,4] = 3
    x[0,5] = 4
    x[0,6] = 5
    x[0,7] = 7
    x[0,8] = 10
    x[0,9] = 30

    x[1,0] = 1
    x[1,1] = 0.9993
    x[1,2] = 0.995
    x[1,3] = 0.982
    x[1,4] = 0.963
    x[1,5] = 0.942
    x[1,6] = 0.9196
    x[1,7] = 0.874
    x[1,8] = 0.806
    x[1,9] = 0.405 

    return x
   
    
def curveInterp(t, curve):
    #zero rates
    tempCurve = -numpy.log(curve[1])/curve[0]
    tempCurve[0] = 0.0000

    if t < 0:
        return 1
    elif t > 30:
        f = numpy.exp(numpy.interp(t, curve[0], numpy.log(curve[1]), right = numpy.log(curve[1,9])*t/30))
        return f
    else:
        f = interpolate.pchip_interpolate(curve[0], tempCurve, t)

    tnew = numpy.exp(-f * t)    
    return tnew

#assume 3m    
def fwdRate(start, curve):
    return (curveInterp(start, curve) / curveInterp(start + 0.25, curve) - 1) / 0.25

def swapRate(start, tenor, curve):
    fix_per = int(tenor * 2)
    flt_per = int(tenor * 4)
    
    flt = 0
    fix = 0
    
    for i in xrange(fix_per):    
        fix = fix + 0.5 * curveInterp(start + i * 0.5, curve)
    
    for j in xrange(flt_per):
        flt = flt + 0.25 * curveInterp(start + j * 0.25, curve) * fwdRate(start + j * 0.25, curve)
    
    return flt / fix * 10000
    
def plotFwd(start, tenor, curve):
    numPer = int(tenor * 4)
    
    fwds = numpy.zeros(numPer)
    dates = numpy.zeros(numPer)
    
    for i in xrange(numPer):
        dates[i] = start + 0.25 * i 
        fwds[i] = fwdRate(dates[i], curve)
    
    pyplot.plot(dates, fwds)
    return
    
def plotZero(start, tenor, curve):
    numPer = int(tenor * 4)
    
    zeros = numpy.zeros(numPer)
    dates = numpy.zeros(numPer)
    
    for i in xrange(numPer):
        dates[i] = start + 0.25 * i 
        zeros[i] = -numpy.log(curveInterp(dates[i], curve))/dates[i]
    
    plot(dates, zeros)
    return

#curveJan = curveBuild('2015-01-30')
#curveJ = yieldCurve('2015-01-30')