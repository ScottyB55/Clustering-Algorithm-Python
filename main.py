# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:26:49 2020

Script analyzes the relationship between stock parameters and overnight returns

Stock parameters determine a stock's location on a cluster's weight curve. The weight curve is defined
by a set of adjustable centers and spreads
A combination of high average returns, low variation, and high number of points across a high number of days
contributes to a cluster's score, which is related to the average daily return we'd expect from buying 
in that cluster

We use log scale to compound through addition, to create more consistent datapoint density, and 
to accurately measure the risk assocatiated with poor returns. Curve seems to yield more consistent
optimization than box


@author: twbur
@author: srburnett
@author rdhurlbu

"""



#                     time            mkt cap     intraday dump         volume      intraday vol     prev night change        day change         vol trend
adjust_cen_spr = [[False, False], [True, True],    [True, True],     [True, True],   [True, True],    [True, True],       [True, True],         [True, True]]

# Steep curve mode, adj_ratio tuned from 1.2 down to 1.05, strat tradeoff 3

"""Stocks that it bought that were actually good:
       CEMI, 314
       AVCT, 317, decent
       SNOA, 313
       RKDA, 312, really good!
       
   Wasteful stocks that it bought:
       HYB , 315
       LATN, 320
       GXGX, 320
   For liquidity purposes, we need to analyze the 15 minute volume data, and avoid
   stocks that have windows of time of little to no volume!
   
   It rarely finds 10 stocks to buy, or invests all the money,
   typically it puts a little bit of money into a crappy stock or two
"""



# Extended out the exponential volume decay a bit... I think it helped!
#               time       mkt cap          intraday dump        volume          intraday vol        prev night change    vol fract         vol trend
cen_spr =       [[15, 50], [16.00, 3.000],  [-0.08, 0.06],       [17.2, 3.52],   [0.07, 0.06],       [0.03, 0.06],        [-3.000, 2.500],  [2.200, 2.000]]

# Optimized up to day 12!
#               [[12, 50], [17.44, 3.0],    [-0.08694, 0.03472], [11.12, 6.083], [0.05139, 0.03472], [0.01272, 0.0864],   [-1.28, 6.221],   [3.156, 1.157]]

cen_spr =       [[15, 50], [16.00, 3.000],  [0.00, 0.06],       [17.2, 3.52],   [0.07, 0.06],       [0.03, 0.06],        [-3.000, 2.500],  [2.200, 2.000]]

#               time       mkt cap          intraday rebound      volume        intraday volit    prev night change     sketch vol fract   vol trend
cen_spr =       [[12, 50], [15.4, 3.0],     [-0.006944, 0.01674], [12.3, 3.52],  [0.128, 0.05],    [-0.01328, 0.03472],  [-0.5, 1.447],     [1.737, 1.157]]

cen_spr =       [[12, 50], [16.0, 3.0],     [-0.006944, 0.01163], [11.6, 2.933], [0.118, 0.05],    [-0.006336, 0.03472], [-0.5, 1.736],     [1.737, 1.157]]
# int gain instead of dump tried, not much difference, probs a tad worse
# Same with close to close tried instead of int dump


#               time       mkt cap          intraday rebound      volume               intraday volit     prev night change     high to open      vol trend
cen_spr =       [[12, 50], [16.0, 3.0],     [-0.006944, 0.01163], [11.6, 2.933],       [0.118, 0.05],     [-0.006336, 0.03472], [0.00, 0.01],     [1.737, 1.157]]

cen_spr =       [[12, 50], [17.73, 4.32],   [-0.004618, 0.01163], [13.2, 6.082],       [0.108, 0.03472],  [-0.007725, 0.04166], [0, 0.005787],    [1.737, 1.157]]

cen_spr =       [[12, 50], [17.73, 4.32],   [0, 0.006944],        [13.2, 10.51],       [0.108, 0.03472],  [0.000607, 0.04166],  [0, 0.006944],    [1.737, 1.157]]


cen_spr =       [[12, 50], [17.73, 4.32],   [0, 0.006944],        [0.00, 0.05],        [0.108, 0.03472],  [0.000607, 0.04166],  [0, 0.006944],    [1.737, 1.157]]


#               time         mkt cap          intraday rebound        prev day change      intraday volit      prev night change     high to open       vol trend
cen_spr =       [[12, 50],   [18.59, 4.32],   [0, 0.006944],          [0.006944, 0.03472], [0.1149, 0.03472],  [0.008939, 0.04999],  [0, 0.006944],     [1.737, 1.157]]

cen_spr =       [[12, 50],   [18.59, 5.184],  [0.001389, 0.006944],   [0.006944, 0.02893], [0.1149, 0.03472],  [0.008939, 0.04999],  [0, 0.006944],     [1.737, 1.157]]


# Try the square root per number day rather than log... log seems better than square root,
# but we could try another form or scaling of log!
# The ln volume trend still seems a little bit better


""" Setup the data source
"""

from _DataSources._DataSourceVol4 import FileMama, DataPointSource

# We start pulling the pickle files from day 0
FileMama = FileMama(0)

# Create a new datasource from the pulled pickles
DataPointSource = DataPointSource(FileMama)

# Allocate memory for datapoints starting 10 (for calculations such as vol that use previous days)
DataPointSource.Allocate_DataPoints(1)

# Fill out the datapoints, with the derived data we actually want
DataPointSource.Update_DataPoints()



""" Setup the cluster analyzer and cache
"""

from _ClusterAnalyzers._CachedClusterAnalyzer2 import CachedClusterAnalyzer

# This determines the tradeoff between the number of datapoints, and the average return confidence
strat_tradeoff = 12    # Squared forced component    180 * 5#3#cen_spr[0][1] * 5#3#2   #-180 * 5
point_n_confidence = 12# Normal Component

# Try making the form e^(-1/x^2)

# This makes it much better! The point_n_confidence kicks in if the strat tradeoff is small, as
# in if we were analyzing a small timeframe, to guarantee that the results aren't random

# The strat tradeoff 

cachedClusterAnalyzer = CachedClusterAnalyzer(DataPointSource, 
                                              cen_spr, 
                                              day_range = (0,9999),
                                              strat_tradeoff = strat_tradeoff, 
                                              point_n_confidence = point_n_confidence,
                                              cluster_type = "steep curve", 
                                              ln_daily_n = True)


""" Let's print out some information about the starting point.
"""

print()
print("Cluster center and spread definition")
DataPointSource.printCenSprArray(cen_spr)
DataPointSource.printCenSpr(cen_spr)
print()
print("Cluster linear range definition")
DataPointSource.printCenSprArray(DataPointSource.convertToLinRange(cen_spr, cachedClusterAnalyzer.dist_2_cutoff * (2 - 0.2) / 2))
print()
print("Cluster statistics")
cachedClusterAnalyzer.Print_Stats(cachedClusterAnalyzer.Score_Curve(cen_spr))
print()



""" Setup the cluster optimizer (tunes centers and spreads towards a local max)
"""

from _ClusterOptimizers._ClusterOptimizer import ClusterOptimizer

""" How much the optimizer scales, in terms of the spread, when stepping through values we are trying out
"""
adj_ratio = 1.2
# Fine tuning
#adj_ratio = 1.1
#adj_ratio = 1.05

clusterOptimizer = ClusterOptimizer(cachedClusterAnalyzer, cen_spr, adjust_cen_spr, adj_ratio)

""" Direct call to cluster optimizer to adjust centers and spreads towards the local max
        clusterOptimizer.Optimize_Step()         # One step across all dimensions
        clusterOptimizer.Optimize_Curve(cen_spr) # Keep stepping until we found the max
"""


from _Strategies._Strategy import Strategy

strategy = Strategy(cachedClusterAnalyzer, buy_threshold = 0.1, buy_limit = 10, max_fraction = 0.33, market_fraction = -9)


from _BackTesters._BackTester import BackTester

""" Backtester can be configured to periodically optimize centers and spreads, and
    track preformance across time
    
    Initialize and run the backtester
"""



backTester =  BackTester(clusterOptimizer, 
                         strategy,
                         start_day = 12,#0#328,250
                         skip_initial_opt = False,#True
                         reopt_period = 99999,
                         start_money = 1000)

# Let her rip across all time
cen_spr[0][1] = 9999

#[[202, 180], [22.72, 12.12], [-0.1041, 0.02757], [8.144, 2.52], [0.06524, 0.03472], [0.03251, 0.1225], [0.07244, 0.0743], [2.512, 2.416]]


# strategy = Strategy(cachedClusterAnalyzer, buy_threshold = 0.1, buy_limit = 10, max_fraction = 0.33, market_fraction = -9)

# backTester =  BackTester(clusterOptimizer, strategy, start_day = 297, start_money = 1000)

# backTester =  BackTester(clusterOptimizer, start_day = 297, start_money = 1000, buy_threshold = 0.1, buy_limit = 10)

# backTester =  BackTester(clusterOptimizer, start_day = 270, start_money = 1000, buy_threshold = 0.1, buy_limit = 10)

# Fixed:
# For some reason, it looks like instantiating another backtester after the run
# doesn't work if we optimized it at the beginning...

# The reason it doesn't work after the run is because the start date gets shifted! It only gets
# points if you go back in time with the start date!

# I think the reason is because it updates the parent cache up until that date so it doesn't cheat!
# ^^^^^^facts, call Update_Parent_Cache before Update_Cache!

# It is doing a bit more than necessary now, its kinda excessive