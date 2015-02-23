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
        self.rnd = Random.randMC(False, False)    

    def spotFwd(self,t):
        lnP = lambda x: numpy.log(yieldCurve.curveInterp(x, self.curve))
        return -scipy.misc.derivative(lnP, t)
 
    def G(self,t, T):
        return (1 - numpy.exp(-self.meanSpeed * (T-t))) / self.meanSpeed
        
    def bondPrice(self, x, t, T):
        Pt = yieldCurve.curveInterp(t, self.curve)
        PT = yieldCurve.curveInterp(T, self.curve)        
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

        
    def eulerPath(self, evoTime, strike, paths, steps):
        dt = evoTime / steps
        vol = self.sigma * numpy.sqrt(dt)
        temp = numpy.zeros(steps+1)        
        discount = 0 
        r_last = self.r0        
        
        rVector = numpy.zeros(steps+1)

        answer = 0
        answerVec = numpy.zeros(steps+1)


        answerVec[0] = self.r0

        #build r0 vector
        for i in xrange(steps+1):
            t = (i)*dt
            rVector[i] = self.spotFwd(i*dt)
            temp[i] = self.sigma**2 / (2 * self.meanSpeed) * (1 - numpy.exp(-2 * self.meanSpeed * t))

        for j in xrange(paths):
            discount = 0         
            r_last = self.r0
            rndNumbers = self.rnd.genNormal(steps)
            x_last = 0                
                
            for i in xrange(steps):
                x_next = x_last +temp[i+1]-temp[i] - self.meanSpeed * x_last * dt + self.sigma * numpy.sqrt(dt) * rndNumbers[i]            
                x_last = x_next                
                r_next = x_next + rVector[i+1]
    
                discount += (r_next + r_last)/2
                r_last = r_next
                answerVec[i+1] = r_next
            
            answer += numpy.exp(-discount*dt)
            #self.plotTCurve(r_next, 10, 30)
        
        return answer / paths

    #risk neutral numeraire
    def exactPath(self, evoTime, strike, paths, steps):
        dt = evoTime / steps
        vol =  numpy.sqrt(self.sigma**2 / (2 * self.meanSpeed) * (1 - numpy.exp(-2 * self.meanSpeed * dt)))
        #term = self.sigma**2/(3*self.meanSpeed) * (1 - numpy.exp(-3*self.meanSpeed*dt))
        y = numpy.zeros(steps+1)
        yIntegral = numpy.zeros(steps+2)      
        yDoubleIntegral = numpy.zeros(steps+2)
        
        discount = 0 
        r_last = self.r0        
        
        rVector = numpy.zeros(steps+1)

        answer = 0
        

        #precompute step
        for i in xrange(steps+1):
            t = (i+1) * dt            
            s = (i) * dt
            rVector[i] = self.spotFwd(i*dt)
            y[i] = self.sigma**2/(2 * self.meanSpeed) * (1 - numpy.exp(-2*self.meanSpeed*s))
            yIntegral[i] = self.sigma**2/(2*self.meanSpeed**2) * ((1 - numpy.exp(-self.meanSpeed*dt)) + (numpy.exp(-2*self.meanSpeed*t) - numpy.exp(-self.meanSpeed*t - self.meanSpeed*s)))
            
            yTerm1 = t - s
            yTerm2 = 1 / self.meanSpeed * (numpy.exp(-self.meanSpeed * (t-s)) - 1)
            yTerm3 =-1 / (2 * self.meanSpeed) * (numpy.exp(-2*self.meanSpeed*t) - numpy.exp(-2 * self.meanSpeed * s))
            yTerm4 = 1 / (self.meanSpeed) * (numpy.exp(-self.meanSpeed*t - self.meanSpeed*s) - numpy.exp(-2*self.meanSpeed*s))
            #yDoubleIntegral[i] = (yIntegral[i] + yIntegral[i+1])/2 * dt
            yDoubleIntegral[i] = self.sigma**2 / (2*self.meanSpeed**2) * (yTerm1 + yTerm2 + yTerm3 + yTerm4)

        #for i in xrange(steps):        
            #print yIntegral[i], (yIntegral[i] + yIntegral[i+1])/2 * dt, yDoubleIntegral[i]

        for j in xrange(paths):
            discount = 0         
            r_last = self.r0
            rndNumbers = self.rnd.genNormal(steps)
            rndNumbers2 = self.rnd.genNormal(steps)
            x_last = 0
            I_last = 0
            
            
            for i in xrange(steps):
                x_next = x_last * numpy.exp(-self.meanSpeed * dt) + yIntegral[i] + vol * rndNumbers[i]

                IVariance = 2 * yDoubleIntegral[i] - y[i] * self.G(i*dt, (i+1)*dt)**2
                I_next = I_last - x_last * self.G(i*dt, (i+1)*dt) - yDoubleIntegral[i] + numpy.sqrt(IVariance) * rndNumbers2[i]
                #print I_next, I_last, 2 * yDoubleIntegral[i], y[i] * self.G(i*dt, (i+1)*dt)**2
                x_last = x_next 
                I_last = I_next
                
                r_next = x_next + rVector[i+1]
                #term = r_last * numpy.exp(-self.meanSpeed * dt) + rVector[i+1] - rVector[i] * numpy.exp(-self.meanSpeed*dt)
                        
                #r_next = term + vol * rndNumbers[i]
                
                #print r_last, term, vol, rVector[i], rVector[i+1], r_next
                discount += (r_next + r_last)/2
                r_last = r_next
            
            #answer += numpy.exp(-discount*dt)
            answer += numpy.exp(I_next)
        
            #self.plotTCurve(r_next, 10, 30)
        
        return answer / paths

model = hullWhite(.1, 0.01, curveJan)

#spot curve
#model.plotTCurve(model.spotR(0), 0, 30)
endDate = 5
print "Analytic = ", yieldCurve.curveInterp(endDate, curveJan)
print "Implied =  ", model.bondPrice(0, 0, endDate)
print "Euler =    ", model.eulerPath(endDate, 0, 5, 200)
print "Exact =    ", model.exactPath(endDate, 0, 100, 1)
