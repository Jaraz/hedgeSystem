# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:34:12 2015

@author: jaycw_000
"""

from __future__ import division
import math
import Random
import numpy
import scipy
from scipy.stats import ncx2
from scipy.special import ive
from matplotlib import pyplot

def logiv(v, z):
    return numpy.log(ive(v,z)) + z

def ncxlogpdf(x, k, l):
    return - math.log(2.) -(x+l)/2. + (k/4. -1./2) * (numpy.log(x) - numpy.log(l)) + logiv(k/2-1, (l * x)**.5)

def ncxpdf(x, k, l):
    return numpy.exp(ncxlogpdf(x, k, l))

pyplot.rcParams.update({'font.size': 26})

def bsOption(expiry, strike, fwd, vol, optType):
    d1 = (numpy.log(fwd / strike) + 0.5 * vol * vol * expiry) / (vol * numpy.sqrt(expiry))
    d2 = d1 - vol * numpy.sqrt(expiry)

    if optType == "call":
        answer = fwd * Random.vecCDF(d1) - strike * Random.vecCDF(d2)
    if optType == "put":
        answer = strike * Random.vecCDF(-d2) - fwd * Random.vecCDF(-d1)
    
    return answer
    
def bsShiftOption(expiry, strike, fwd, vol, shift, optType):
    return shiftLNOption(expiry, strike, fwd, vol, shift * vol, optType)

def normOption(expiry, strike, fwd, vol, optType):
    d1 = (fwd - strike) / (vol * numpy.sqrt(expiry))
    d2 = -(fwd - strike) / (vol * numpy.sqrt(expiry))
    nd1Prime = 1/numpy.sqrt(2 * numpy.pi) * numpy.exp(-d1**2/2)
    nd2Prime = 1/numpy.sqrt(2 * numpy.pi) * numpy.exp(-d2**2/2)
    
    if optType == "call":
        answer = vol * numpy.sqrt(expiry) * (d1 * Random.vecCDF(d1) + nd1Prime)
    if optType == "put":
        answer = vol * numpy.sqrt(expiry) * (d2 * Random.vecCDF(d2) + nd2Prime)
    
    return answer

def shiftLNOption(expiry, strike, fwd, vol1, vol2, optType):
    d1 = (numpy.log((vol1 * fwd + vol2)/(vol1 * strike + vol2)) + 0.5 * vol1 * vol1 * expiry) / (vol1 * numpy.sqrt(expiry))
    d2 = d1 - vol1 * numpy.sqrt(expiry)
    
    if optType == "call":
        answer = (fwd + vol2 / vol1) * Random.vecCDF(d1) - (strike + vol2/vol1) * Random.vecCDF(d2)
    if optType == "put":
        answer = -(fwd + vol2 / vol1) * Random.vecCDF(-d1) + (strike + vol2/vol1) * Random.vecCDF(-d2)

    return answer

def normVol(expiry, strike, fwd, optType, optPrice):
    optFunc = lambda x: normOption(expiry, strike, fwd, x, optType) - optPrice
    return scipy.optimize.brentq(optFunc, a= 0.0001, b=10000)

def cevIntegrate(expiry, strike, fwd, vol, beta, optType):
    nu = 0.5 / (1.0 - beta)
    term1 = (4 * nu * fwd) / (vol * vol * expiry)
    term2 = (4 * nu * nu) / (vol * vol * expiry)
    term3 = 2 * nu + 2
    term4 = (4 * nu * nu * numpy.power(fwd, 1.0 / nu)) / (vol * vol * expiry)

    if optType == "call":
        integrand = lambda f: (f - strike) * term1 * numpy.power(f,1.0/nu-2) * ncxpdf(term2 * numpy.power(f,1/nu),term3,term4)

        answer = scipy.integrate.quad(integrand, strike, numpy.inf, epsabs=1.19e-19, epsrel=1.19e-19, limit=100)[0]
    
    if optType == "put":
        integrand = lambda f: (strike - f) * term1 * numpy.power(f,1.0/nu-2) * ncxpdf(term2 * numpy.power(f,1/nu),term3,term4)
        answer =   (scipy.integrate.quad(integrand, 0, strike, epsabs=1.19e-19, epsrel=1.19e-19, limit=100)[0] + 
                    max(strike,0) * 1/scipy.special.gamma(nu) * scipy.special.gammaincc(nu, (2 * nu**2 * numpy.power(fwd,1/nu)) / (vol**2 * expiry)))
    
    return answer

def cevOption(expiry, strike, fwd, vol, beta, optType):
    if beta == 1:
        return bsOption(expiry, strike, fwd, vol, optType)
    
    nu = 0.5 / (1.0 - beta)
    
    densityFunc = lambda x, r, lmbda: ( 0.5 * numpy.power( (x/lmbda), ((r-2)*0.25) ) * 
                                        numpy.exp(-0.5*(x+lmbda)) * 
                                        scipy.special.iv(0.5*(r-2),numpy.sqrt(lmbda*x)))
    
    term11 = (4 * numpy.power(nu,2) * numpy.power(strike,1.0/nu)) / (numpy.power(vol,2) * expiry)
    term12 = 2 * nu + 2
    term13 = (4 * numpy.power(nu,2) * numpy.power(fwd,(1.0/nu))) / (numpy.power(vol,2) * expiry)
    term22 = 2 * nu
    
    #int1 = scipy.integrate.quad(densityFunc, 0, term11, args=(term12, term13))[0]
    #int2 = scipy.integrate.quad(densityFunc, 0, term13, args=(term22, term11))[0]
    
    #testX = 816
    #print ( densityFunc(testX, term12, term13), 
    #        0.5 * numpy.power( (testX/term13), ((term12-2)*0.25) ),
    #        numpy.exp(-0.5*(testX+term13)),
    #        scipy.special.iv(0.5*(term12-2),numpy.sqrt(term13*testX)) )
    int1 = ncx2.cdf(term11, term12, term13)
    int2 = ncx2.cdf(term13, term22, term11)
    print fwd, strike, int1, int2
    
    return fwd * (1 - int1) - strike * int2
    

def optionGraph(expiry, strikeVec, fwd, vol, shift1, shift2, optType, graphType):
    #bs     = bsOption(expiry, strikeVec, fwd, vol1, optType)
    atmNormPrice = normOption(expiry, fwd, fwd, vol, optType)
    norm   = normOption(expiry, strikeVec, fwd, vol, optType)
    shiftVol1 = volShiftConvert(expiry, fwd, fwd, shift1, optType, atmNormPrice)
    shiftVol2 = volShiftConvert(expiry, fwd, fwd, shift2, optType, atmNormPrice)

    shiftAns1  = bsShiftOption(expiry, strikeVec, fwd, shiftVol1, shift1, optType)
    shiftAns2  = bsShiftOption(expiry, strikeVec, fwd, shiftVol2, shift2, optType)

    #pyplot.plot(strikeVec, bs)
    if graphType == "price":    
        pyplot.plot(strikeVec, norm)
        pyplot.plot(strikeVec, shiftAns1)  
        pyplot.plot(strikeVec, shiftAns2)
    
    if graphType == "vol":
        pyplot.plot(strikeVec, numpy.repeat(vol, strikeVec.size), label="normATM")
        pyplot.plot(strikeVec, normConvert(expiry, strikeVec, fwd, optType, shiftAns1), label="shift1")
        pyplot.plot(strikeVec, normConvert(expiry, strikeVec, fwd, optType, shiftAns2), label="shift2")
        pyplot.legend()

def normConvert(expiry, strike, fwd, optType, optPrice):
    answer = numpy.zeros(strike.size)    
    for i in range(1,strike.size):
        answer[i] = normVol(expiry, strike[i], fwd, optType, optPrice[i])
    
    return answer

#find shift vol for given zeroShift level to match given option price
def volShiftConvert(expiry, strike, fwd, zeroShift, optType, optPrice):
    optFunc = lambda x: bsShiftOption(expiry, strike, fwd, x, zeroShift, optType) - optPrice
    return scipy.optimize.brentq(optFunc, a = 0.0001, b = 100)
    
testVec = numpy.linspace(0,30,41)
shiftVec = numpy.linspace(10, 100, 6)

print cevIntegrate(1, 10, 100, 4, .5, "put")
print normOption(1, 10, 100, 10, "put")
print bsOption(1, 10, 100, .1, "put")
optionGraph(1, testVec, 20, 20, 10, 100, "put", "vol")