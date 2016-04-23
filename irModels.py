# -*- coding: utf-8 -*-
"""
Created on Sun Feb 15 21:50:49 2015

@author: jaycw_000
"""

from __future__ import division
import numpy
import scipy
import yieldCurve
import Random
import matplotlib
from matplotlib import pyplot

class irModel:
    def path(self):
        return 1

class hoLee:
    def __init__(self, sigma, curve):
        self.sigma = sigma
        self.curve = curve
        self.r0 = self.f(0)
        self.rnd = Random.randMC(False, False)
        
    def updateR(self, r):
        self.r0 = r

    def updateTheta(self, sigma):
        self.sigma = sigma
        
    def f(self, t):
        lnP = lambda x: numpy.log(yieldCurve.curveInterp(x, self.curve))
        return -scipy.misc.derivative(lnP, t)
    
    def f2(self, t):
        return scipy.misc.derivative(self.f, t)

    def bondPrice(self, r, t, T):
        pt = yieldCurve.curveInterp(t, self.curve)
        pT = yieldCurve.curveInterp(T, self.curve)
        
        term1 = -(r - self.f(t)) * (T - t)
        term2 = 0.5 * self.sigma**2 * t * (T-t)**2
        
        return pT / pt * numpy.exp(term1 - term2)

    def plotTCurve(self, r, t, T):
        numFwds = int((T-t)*4 - 1)
        dates = numpy.zeros(numFwds+1)        
        fwds = numpy.zeros(numFwds)

        dates[0] = t

        for i in xrange(numFwds):
            dates[i+1] = t + (i+1) * 0.25
            fwds[i] = 1 / 0.25 * (self.bondPrice(r, t, dates[i]) / self.bondPrice(r, t, dates[i+1]) - 1)
        
        plot(dates[0:numFwds],fwds)
    
    def plotSpotCurve(self):
        numPer = 120
        
        fwds = numpy.zeros(numPer-1)
        dates = numpy.zeros(numPer)
        bonds = numpy.zeros(numPer)
        
        start = 0        

        #generate bonds        
        for i in xrange(numPer):
            dates[i] = start + 0.25 * i 
            bonds[i] = numpy.exp(-self.r0 * dates[i] + 1/6 * self.sigma**2 * dates[i]**3) 
        
        for j in xrange(numPer-1):        
            fwds[j] = (bonds[j] / bonds[j+1] - 1) / 0.25
        
        plot(dates[0:119], fwds)
    
    def path(self, evoTime, steps):
        dt = evoTime / steps
        rndNumbers = self.rnd.genNormal(steps)
        r_last = self.r0        
        
        bond_price = 0        
        discount = 0        
        
        for i in xrange(steps):
            t = dt * i

            r_next = r_last + self.f(t+dt) - self.f(t) + self.sigma**2 / 2 * ((t+dt)**2 - t**2) + self.sigma*sqrt(dt)*rndNumbers[i]
            discount += r_last
            r_last = r_next

        #self.plotTCurve(r_next, t, 30)
        #return exp(-discount*dt) * 1/0.25 * (1/self.bondPrice(r_next, 5, 5.25)-1)
        return exp(-discount*dt) * self.bondPrice(r_next,5,5)
    
    
    
##Regular Vasicek, not extended
class vasicek:
    def __init__(self, r0, meanLevel, meanSpeed, sigma, curve):
        self.sigma = sigma
        self.curve = curve
        self.r0 = r0
        self.meanLevel = meanLevel
        self.meanSpeed = meanSpeed
        self.rnd = Random.randMC(False, False)
    
    def bondPrice(self, r, t, T):
        B = (1 - numpy.exp(-self.meanSpeed * (T-t))) / self.meanSpeed
        A = (self.meanLevel - self.sigma**2/(2*self.meanSpeed**2)) * (B - (T-t)) - (self.sigma**2 * B**2)/(4*self.meanSpeed)
        
        P = numpy.exp(A - B * r)
        
        return P
    
    def optionPrice(self, t, T, strike, S):
        B = (1 - numpy.exp(-self.meanSpeed * (T-t))) / self.meanSpeed
        A = (self.meanLevel - self.sigma**2/(2*self.meanSpeed**2)) * (B - (T-t)) - (self.sigma**2 * B**2)/(4*self.meanSpeed)
        Bts = (1 - numpy.exp(-self.meanSpeed * (S-T))) / self.meanSpeed        
        strikeConvert = 1 / (1 + strike * 0.25)
        
        sigmaP = self.sigma * numpy.sqrt((1 - numpy.exp(-2*self.meanSpeed*(T-t))) / (2 * self.meanSpeed)) * Bts
        h = 1 / sigmaP * numpy.log(self.bondPrice(self.r0, t, S) / (self.bondPrice(self.r0, t, T)*strikeConvert)) + sigmaP / 2

        return self.bondPrice(self.r0, t, S) * scipy.stats.norm.cdf(h) - strikeConvert * self.bondPrice(self.r0, t, T) * scipy.stats.norm.cdf(h-sigmaP)
        
    
    def plotSpotCurve(self):
        numPer = 120
        
        fwds = numpy.zeros(numPer-1)
        dates = numpy.zeros(numPer)
        bonds = numpy.zeros(numPer)
        
        start = 0        

        #generate bonds        
        for i in xrange(numPer):
            dates[i] = start + 0.25 * i 
            bonds[i] = self.bondPrice(self.r0, 0, dates[i])
        
        for j in xrange(numPer-1):        
            fwds[j] = (bonds[j] / bonds[j+1] - 1) / 0.25
        
        plot(dates[0:119], fwds)
        
    def path(self, evoTime, strike, steps):
        dt = evoTime / steps
        rndNumbers = self.rnd.genNormal(steps)
        r_last = self.r0        
        
        bond_price = 0        
        discount = 0        
        
        for i in xrange(steps):
            t = dt * i

            r_next = r_last + self.meanSpeed * (self.meanLevel - r_last) * dt + self.sigma*numpy.sqrt(dt)*rndNumbers[i]
            discount += r_last
            r_last = r_next
        
        #self.plotTCurve(r_next,1,30)
        #return exp(-discount*dt) * self.bondPrice(r_next,1,1)
        strikeConvert = 1 / (1 + strike * 0.25)
        return exp(-discount*dt) * max(self.bondPrice(r_next, t, t+0.25) - strikeConvert, 0)

    def plotTCurve(self, r, t, T):
        numFwds = int((T-t)*4 - 1)
        dates = numpy.zeros(numFwds+1)        
        fwds = numpy.zeros(numFwds)

        dates[0] = t

        for i in xrange(numFwds):
            dates[i+1] = t + (i+1) * 0.25
            fwds[i] = 1 / 0.25 * (self.bondPrice(r, t, dates[i]) / self.bondPrice(r, t, dates[i+1]) - 1)
        
        plot(dates[0:numFwds],fwds)


#Implemented from Andersen and Piterbarg
class hullWhite:
    def __init__(self, meanSpeed, sigma, curve):
        self.sigma = sigma
        self.curve = curve
        self.r0 = self.spotFwd(0)
        self.meanSpeed = meanSpeed
        self.rnd = Random.randMC(True, False)    

    def spotFwd(self,t):
        lnP = lambda x: numpy.log(self.curve.discFact(x))
        return -scipy.misc.derivative(lnP, t)
 
    def G(self,t, T):
        return (1 - numpy.exp(-self.meanSpeed * (T-t))) / self.meanSpeed
    
    def y(self, t):
        return self.sigma**2 / (2 * self.meanSpeed) * (1 - numpy.exp(-2 * self.meanSpeed * t))
    
    def IVar(self, ti, ti1):
        a = self.meanSpeed
        
        term1 = ti1 - ti
        term2 = 1 / a     * (numpy.exp(-a*(ti1-ti)) - 1)
        term3 = 1 / (2*a) * (numpy.exp(-2*a*ti1) - numpy.exp(-2*a*ti))
        term4 = 1 / a     * (numpy.exp(-a*(ti1+ti)) - numpy.exp(-2*a*ti))

        return self.sigma**2/(a**2) * (term1 + term2 - term3 + term4) - self.y(ti) * self.G(ti, ti1)**2
    
    def bondPrice(self, x, t, T):
        Pt = self.curve.discFact(t)
        PT = self.curve.discFact(T)
        GVar = self.G(t,T)
        y = self.sigma**2 / (2 * self.meanSpeed) * (1 - numpy.exp(-2*self.meanSpeed*t))

        P = PT/Pt * numpy.exp(-x * GVar - 0.5 * y * GVar**2)
        
        return P

    def plotTCurve(self, x, t, T):
        numFwds = int((T-t)*4 - 1)
        dates = numpy.zeros(numFwds+1)        
        fwds = numpy.zeros(numFwds)

        dates[0] = t

        for i in xrange(numFwds):
            dates[i+1] = t + (i+1) * 0.25
            fwds[i] = 1 / 0.25 * (self.bondPrice(x, t, dates[i]) / self.bondPrice(x, t, dates[i+1]) - 1)
        
        pyplot.plot(dates[0:numFwds],fwds)

    def cov(self, ti, ti1, sigma, a):
        term1 = 1 - numpy.exp(-a * (ti1-ti))
        term2 = numpy.exp(-2*a*(ti1-ti)) - numpy.exp(-a*(ti1-ti))
    
        return -sigma**2 / (2*a*a)*(term1+term2)
        
    def capFloorPricer(self, expiry, strike, optType):
        #convert strike from bps to DF
        adjStrike = 1 / (1 + 0.25 * strike)
        
        ptMat = self.curve.discFact(expiry + 0.25)
        pt    = self.curve.discFact(expiry)
        
        h = (1 - numpy.exp(-self.meanSpeed*0.25)) / self.meanSpeed
        
        vol = self.sigma * h * numpy.sqrt(1 / (2*self.meanSpeed) * (1 - numpy.exp(-2*self.meanSpeed*expiry)))
        
        d1 = 1 / vol * numpy.log(ptMat / (pt * adjStrike)) + vol / 2
        d2 = d1 - vol
        
        if optType == "Cap":
            return ptMat * Random.normCDF(d1) - adjStrike * pt * Random.normCDF(d2)
        elif optType == "Floor":
            return adjStrike * pt * Random.normCDF(-d2) - ptMat * Random.normCDF(-d1)
        
        return 999

    def swapPricer(self, x, start, tenor, strike, payRec, debug = False):
        #gen dates
        fix = 0
        flt = 0
        
        for i in xrange(1,tenor*2+1):
            temp = start + 0.5 * i
            fix += 0.5 * self.bondPrice(x, start, temp) * strike
            if debug == True:
                print temp, fix
                
                    
        flt += 1 - self.bondPrice(x, start, start+tenor)
        
        if payRec == "Pay":
            return fix - flt
        else:
            return flt - fix
        
    def eulerPath(self, evoTime, strike, paths, steps):
        dt = evoTime / steps
        vol = self.sigma * numpy.sqrt(dt)
        y = numpy.zeros(steps+1)        
        rndMatrix = self.rnd.genNormalMatrix(paths*steps, paths, steps)
        xAnswer = numpy.zeros(paths)        
        xLast   = numpy.zeros(paths)
        iAnswer = numpy.zeros(paths)
        iLast   = numpy.zeros(paths)        

        answer = 0

        #build r0 vector
        for i in xrange(steps+1):
            t = (i)*dt
            y[i] = self.sigma**2 / (2 * self.meanSpeed) * (1 - numpy.exp(-2 * self.meanSpeed * t))

        for i in xrange(steps):
            xAnswer = numpy.exp(-self.meanSpeed*dt) * xLast + (1 - numpy.exp(-self.meanSpeed*dt))/self.meanSpeed * y[i] + vol * rndMatrix[:,i]            
            iAnswer = iLast - xAnswer * dt

            xLast   = xAnswer
            iLast   = iAnswer

        #answer = numpy.exp(iAnswer) * numpy.maximum(self.bondPrice(xAnswer, evoTime, evoTime+0.25) - 1 / (1 + 0.25 * strike), 0)
        #answer = numpy.exp(iAnswer) * self.swapPricer(xAnswer, 10, 10, strike, "Pay")
        answer = numpy.exp(iAnswer) * self.bondPrice(xAnswer, 10,20)
         
        return self.curve.discFact(evoTime) * answer.mean() * 10000

    #risk neutral numeraire
    def exactPath(self, evoTime, strike, paths, steps):
        dt = evoTime / steps
        IVariance = numpy.zeros(steps)
        IVol = numpy.zeros(steps)
        y = numpy.zeros(steps)
        yIntegral = numpy.zeros(steps)      
        covariance = numpy.zeros(steps)
        corr = numpy.zeros(steps)
        rndCorr = numpy.zeros(steps)
        gVec = numpy.zeros(steps)
        rndMatrix1 = self.rnd.genNormalMatrix(paths*steps, paths, steps)
        rndMatrix2 = self.rnd.genNormalMatrix(paths*steps, paths, steps)
        answer = 0

        expMeanSpeed = numpy.exp(-self.meanSpeed * dt)
        vol =  numpy.sqrt(self.sigma**2 / (2 * self.meanSpeed) * (1 - numpy.exp(-2 * self.meanSpeed * dt)))

        #precompute steps
        for i in xrange(steps):
            t = (i+1) * dt            
            s = (i) * dt
            y[i] = self.y(s)
            yIntegral[i] = self.sigma**2/(2*self.meanSpeed**2) * ((1 - numpy.exp(-self.meanSpeed*dt)) + (numpy.exp(-2*self.meanSpeed*t) - numpy.exp(-self.meanSpeed*t - self.meanSpeed*s)))

            IVariance[i]  = self.IVar(s, t)
            IVol[i] = numpy.sqrt(IVariance[i])
            covariance[i] = self.cov(s, t, self.sigma, self.meanSpeed)
            corr[i] = covariance[i] / (vol * IVol[i])
            gVec[i] = self.G(s, t)

        rndCorr = rndMatrix1 * corr + rndMatrix2 * numpy.sqrt(1-corr**2)

        xAnswer = numpy.zeros(paths)        
        xLast   = numpy.zeros(paths)
        iAnswer = numpy.zeros(paths)
        iLast   = numpy.zeros(paths)        

        for i in xrange(steps):
            xAnswer = xLast * expMeanSpeed + yIntegral[i] + vol * rndMatrix1[:,i]
            iAnswer = iLast - xLast * gVec[i] -(IVariance[i] + y[i] * gVec[i]) / 2 + IVol[i]*rndCorr[:,i]
            
            xLast = xAnswer
            iLast = iAnswer
        answer = numpy.exp(iAnswer)
        #answer = numpy.exp(iAnswer) * numpy.maximum(self.bondPrice(xAnswer, evoTime, evoTime+0.25) - 1 / (1 + 0.25 * strike), 0)
        #answer = numpy.exp(iAnswer) * self.swapPricer(xAnswer, 10, 10, strike, "Pay")
        #answer = numpy.exp(iAnswer) * self.bondPrice(xAnswer, 10,10.25)
        return answer.mean()
        return self.curve.discFact(evoTime) * answer.mean() * 10000

#model = hullWhite(0.02, 0.003, curveJ)

#spot curve
#model.plotTCurve(0, 0, 30)
#endDate = 10
#print "Analytic    =  ", yieldCurve.swapPricer(10,10,0.023096, curveJan) * 10000
#print 'Anal test   =  ', model.swapPricer(0, 10, 10, 0.023096, "Pay", True) * 10000
#print "Analytic    =  ", model.capFloorPricer(endDate, 0.020919083, "Cap") * 10000
#print "Analytic    =  ", model.bondPrice(0,0,10.25) * 10000
#print "Euler MC    =  ", model.eulerPath(endDate, 0.023096, 10000, 100)
#print "Monte Carlo =  ", model.exactPath(endDate, 0.023096, 10000, 1)