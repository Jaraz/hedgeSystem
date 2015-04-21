# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 19:31:45 2015

@author: jaycw_000
"""

from __future__ import division
import numpy


class tree:
    def __init__(self):
        self.recomb = True
    
    #def createTree(self, S, r, vol, T, steps):
        
    
    def bsTree(self, S, r, vol, optionType, steps):
        T = optionType.returnExpiry()
        dt = T / steps
        disc = numpy.exp(-r * dt)
        
        #Build Tree
        tree = [numpy.zeros(1)]
        answer = [numpy.zeros(1)]
        
        for i in xrange(1,steps):
            tree.append(numpy.zeros(i+1))
            answer.append(numpy.zeros(i+1))

        tree[0][0] = numpy.log(S)
        up = (r - 0.5 * vol**2) * dt + vol * numpy.sqrt(dt)
        down = (r - 0.5 * vol**2) * dt - vol * numpy.sqrt(dt)
                
        for i in xrange(1,steps):
            for j in xrange(i,0,-1):
                tree[i][j] = tree[i-1][j-1] + up
            tree[i][0] = tree[i-1][0] + down
        
        #Price option
        answer[steps-1] = optionType.payoff(numpy.exp(tree[steps-1]))
        
        for i in xrange(steps - 2, -1, -1):
            answer[i] = disc * 0.5 *(answer[i+1][1:] + answer[i+1][:-1])
            if optionType.returnEuropean() == False:
                answer[i] = numpy.maximum(optionType.payoff(numpy.exp(tree[i])),answer[i])

        return tree, answer
        
testTree = tree()
putOption = put(100, 1, 1, False)
temp, tempA = testTree.bsTree(100, .1, 0.3, putOption, 50)
print tempA[0][0]