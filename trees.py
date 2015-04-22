# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 19:31:45 2015

@author: jaycw_000
"""

from __future__ import division
import numpy

class treeType:
    def __init__(self):
        self.desc = "tree"

class logTree(treeType):
    def returnUp(self, S, r, vol, dt):
        up = (r - 0.5 * vol**2) * dt + vol * numpy.sqrt(dt)
        return up
    
    def returnDown(self, S, r, vol, dt):
        down = (r - 0.5 * vol**2) * dt - vol * numpy.sqrt(dt)
        return down

class normalTree(treeType):
    def returnUp(self, S, r, vol, dt):
        up = r * dt + numpy.log(vol * numpy.sqrt(dt))
        return up

    def returnDown(self, S, r, vol, dt):
        up = r * dt - numpy.log(vol * numpy.sqrt(dt))
        return up        

class tree:
    def __init__(self, treeModel):
        self.recomb = True
        self.treeModel = treeModel   #Type tree
    
    def createTree(self, S, r, vol, T, steps):
        #Build Log Tree
        dt = T / steps
        tree = [numpy.zeros(1)]        

        for i in xrange(1,steps):
            tree.append(numpy.zeros(i+1))        

        tree[0][0] = numpy.log(S)
        up = self.treeModel.returnUp(S, r, vol, dt)
        down = self.treeModel.returnDown(S, r, vol, dt)
                
        for i in xrange(1,steps):
            tree[i][1:] = tree[i-1] + up
            tree[i][0] = tree[i-1][0] + down
        
        return tree
            
    def priceTree(self, S, r, vol, optionType, steps):
        T = optionType.returnExpiry()
        dt = T / steps
        disc = numpy.exp(-r * dt)
        tree = self.createTree(S, r, vol, T, steps)
        
        #Build answer tree
        answer = [numpy.zeros(1)]
        for i in xrange(1,steps):
            answer.append(numpy.zeros(i+1))

        #Price option - last step first
        answer[steps-1] = optionType.payoff(numpy.exp(tree[steps-1]))
        
        for i in xrange(steps - 2, -1, -1):
            answer[i] = disc * 0.5 *(answer[i+1][1:] + answer[i+1][:-1])
            if optionType.returnEuropean() == False:
                answer[i] = numpy.maximum(optionType.payoff(numpy.exp(tree[i])),answer[i])

        return tree, answer

logType = logTree()
normType = normalTree()
testTree = tree(logType)
testTree2 = tree(normType)
putOption = put(100, 1, 1, False)
temp, tempA = testTree.priceTree(100, 0.05, 0.1, putOption, 200)
temp2, tempA2 = testTree2.priceTree(100, 0.05, 0.1, putOption, 200)
print tempA[0][0], tempA2[0][0]