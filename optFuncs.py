# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 09:35:47 2015

@author: jaycw_000
"""

from __future__ import division
import math
import Random

def forwardVal(T, K, S, r):
    return(S - K*math.exp(-r*T))

def optionVal(T, K, callPut, S, r, vol):
    d1 = (math.log(S / K) + (r + 0.5 * vol * vol) * T) / (vol * T ** 0.5)
    d2 = d1 - vol * T ** 0.5

    if callPut == "call":
        answer = S * Random.normCDF(d1) - K * math.exp(-r*T) * Random.normCDF(d2)
    elif callPut == "put":
        answer = K * math.exp(-r*T) * Random.normCDF(-d2) - S * Random.normCDF(-d1)
    elif callPut == "straddle":
        answer = S * Random.normCDF(d1) - K * math.exp(-r*T) * Random.normCDF(d2) + K * math.exp(-r*T) * Random.normCDF(-d2) - S * Random.normCDF(-d1)      
    elif callPut == "digCall":
        answer = math.exp(-r*T) * Random.normCDF(d2)
    elif callPut == "digPut":
        answer = math.exp(-r*T) * Random.normCDF(-d2)
    elif callPut == "stock":
        answer = S
    elif callPut == "fwd":
        answer = forwardVal(T,K,S,r)
    return answer

def calcGamma(S, r, vol, portfolio):
    answer = 0
    port = portfolio.returnSec()
    
    for i in xrange(len(port)):
        T = port[i].returnExpiry()
        K = port[i].returnStrike()
        
        answer += port[i].returnUnits() * optionGamma(T, K, port[i].returnDesc(), S, r, vol)
    
    return answer

def calcVega(S, r, vol, portfolio):
    answer = 0
    port = portfolio.returnSec()
    
    for i in xrange(len(port)):
        T = port[i].returnExpiry()
        K = port[i].returnStrike()
        
        answer += port[i].returnUnits() * optionVega(T, K, port[i].returnDesc(), S, r, vol)
    
    return answer
    
def calcDelta(S, r, vol, portfolio):
    answer = 0
    port = portfolio.returnSec()
    
    for i in xrange(len(port)):
        T = port[i].returnExpiry()
        K = port[i].returnStrike()
        
        answer += port[i].returnUnits() * optionDelta(T, K, port[i].returnDesc(), S, r, vol)
        
    return answer

def calcPrice(S, r, vol, portfolio):
    answer = 0
    port = portfolio.returnSec()
    
    for i in xrange(len(port)):
        T = port[i].returnExpiry()
        K = port[i].returnStrike()
        
        answer += port[i].returnUnits() * optionVal(T, K, port[i].returnDesc(), S, r, vol)
    
    return answer

def calcPayoff(S, portfolio):
    answer = 0
    port = portfolio.returnSec()
    
    for i in xrange(len(port)):
        answer += port[i].payoff(S)
    return answer

def optionDelta(T, K, callPut, S, r, vol):
    d1 = (math.log(S / K) + (r + 0.5 * vol * vol) * T) / (vol * T ** 0.5)
    d2 = d1 - vol * T ** 0.5
    
    if callPut == "call":
        answer = Random.normCDF(d1)
    if callPut == "put":
        answer = Random.normCDF(d1) - 1
    if callPut == "digCall":
        answer = math.exp(-r * T) * nPrime(d2) / (S * vol * math.sqrt(T))
    if callPut == "digPut":
        answer = - math.exp(-r * T) * nPrime(d2) / (S * vol * math.sqrt(T))
    if callPut == "fwd":
        answer = math.exp(-r * T)
    if callPut == "stock":
        answer = 1
        
    return answer
  
def nPrime(x):
    return math.exp(-x**2 / 2) / math.sqrt(2 * math.pi)
 
def optionGamma(T, K, callPut, S, r, vol):
    d = (math.log(S / K) + (r + 0.5 * vol * vol) * T) / (vol * T ** 0.5)
    d2 = d - vol * T ** 0.5
    #nPrime = math.exp(-d ** 2 / 2) / math.sqrt(2 * pi)
    
    if callPut == "call":
        answer = nPrime(d) / (S * vol * math.sqrt(T))
    elif callPut == "put":
        answer = nPrime(d) / (S * vol * math.sqrt(T))        
    elif callPut == "digCall":
        answer = -math.exp(-r * T) * d * nPrime(d2) / (vol*vol * S*S * T)
    elif callPut == "digPut":
        answer = math.exp(-r * T) * d * nPrime(d2) / (vol**2 * S**2 * T)
    elif callPut == "stock":
        answer = 0
    elif callPut == "fwd":
        answer = 0
        
    return answer

def optionVega(T, K, callPut, S, r, vol):
    d = (math.log(S / K) + (r + 0.5 * vol * vol) * T) / (vol * T ** 0.5)
    d2 = d - vol * math.sqrt(T)
    
    if callPut == "call":
        answer = S * math.sqrt(T) * nPrime(d)
    elif callPut == "put":
        answer = S * math.sqrt(T) * nPrime(d)
    elif callPut == "digCall":
        answer = - math.exp(-r * T) * nPrime(d2) * d / vol
    elif callPut == "digPut":
        answer = - math.exp(-r * T) * nPrime(d2) * d / vol
    elif callPut == "stock":
        answer = 0
    elif callPut == "fwd":
        answer = 0
        
    return(answer)    

def optionTheta(T, K, callPut, S, r, vol):
    d1 = (math.log(S / K) + (r + 0.5 * vol * vol) * T) / (vol * T ** 0.5)
    d2 = d1 - vol * T ** 0.5
    nPrime = math.exp(-d1 ** 2 / 2) / math.sqrt(2 * math.pi)

    if callPut == "call":
        term1 = -r * K * math.exp(-r * T) * Random.normCDF(d2)
        term2 = (vol * S * nPrime) / (2 * math.sqrt(T))
        answer = term1 - term2
    
    if callPut == "put":
        term1 = r * K * math.exp(-r * T) * Random.normCDF(-d2)
        term2 = (vol * S * nPrime) / (2 * math.sqrt(T))
        answer = term1 - term2
        
    return answer
    
def bachPrice(T, K, callPut, S, r, vol):
    d = (S * math.exp(r * T) - K) / (vol * math.sqrt(T))
    
    if callPut == "call":
        answer = (S * math.exp(r * T) - K) * Random.normCDF(d) + vol * math.sqrt(T) * 1 / math.sqrt(2 * math.pi) * math.exp(-d**2/2)
    if callPut == "put":
        answer = 0
    
    return answer