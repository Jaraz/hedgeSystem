# -*- coding: utf-8 -*-
from Quandl import Quandl
from matplotlib import pyplot
import pandas

#quadl authentication token
auth_tok = "EuSxuE2sgxd21EWpW-J9"

class edFuture:
    #load fed funds spot contract
    def __init__(self):
        #create mapping
        self.edDict = {1:"CME_ED1", 
                       2:"CME_ED2", 
                       3:"CME_ED3",
                       4:"CME_ED4",
                       5:"CME_ED5",
                       6:"CME_ED6",
                       7:"CME_ED7",
                       8:"CME_ED8",
                       9:"CME_ED9",
                       10:"CME_ED10",
                       11:"CME_ED11",
                       12:"CME_ED12",
                       13:"CME_ED13",
                       14:"CME_ED14",
                       15:"CME_ED15",
                       16:"CME_ED16"}
        ed1 = Quandl.get("CHRIS/CME_FF1", trim_start = "2005-12-04", authtoken = auth_tok)
        self.df = 100 - pandas.DataFrame(ed1["Settle"])
        self.df.columns = ["CME_FF1"] 
    
    def loadED(self, n):
        label = "CHRIS/" + self.edDict[n]
        self.df[self.edDict[n]] = 100 - Quandl.get(label, trim_start = "2005-12-04", authtoken = auth_tok)["Settle"]
   
    #first 16 contracts
    def loadAll(self):
        for i in range(1,17):
            label = "CHRIS/" + self.edDict[i]
            self.df[self.edDict[i]] = 100 - Quandl.get(label, trim_start = "2005-12-04", authtoken = auth_tok)["Settle"]

        #fix bad data on website
        self.df.ix[1477,15] = (self.df.ix[1476,15] + self.df.ix[1478,15])/2
        self.df.ix[1477,16] = (self.df.ix[1476,16] + self.df.ix[1478,16])/2

        self.df.ix[1475,15] = (self.df.ix[1474,15] + self.df.ix[1476,15])/2
        self.df.ix[1475,16] = (self.df.ix[1474,16] + self.df.ix[1476,16])/2

    #load everything first, x - y
    def spread(self, x,y):
        label  = "spread_" + str(x) + "_" + str(y)
        label1 = self.edDict[y]
        label2 = self.edDict[x]
        
        self.df[label] = self.df[label2] - self.df[label1]

    def runColors(self):
        self.df["whites"] = ( self.df["CME_ED1"] +
                              self.df["CME_ED2"] +
                              self.df["CME_ED3"] +
                              self.df["CME_ED4"]  ) / 4
        
        self.df["reds"]   = ( self.df["CME_ED5"] +
                              self.df["CME_ED6"] +
                              self.df["CME_ED7"] +
                              self.df["CME_ED8"]  ) / 4
        
        self.df["greens"] = ( self.df["CME_ED9"] +
                              self.df["CME_ED10"] +
                              self.df["CME_ED11"] +
                              self.df["CME_ED12"]  ) / 4

        self.df["blues"]  = ( self.df["CME_ED13"] +
                              self.df["CME_ED14"] +
                              self.df["CME_ED15"] +
                              self.df["CME_ED16"]  ) / 4



ed = edFuture()
ed.loadAll()
ed.spread(2,1)
ed.spread(4,1)
ed.spread(8,1)
ed.runColors()