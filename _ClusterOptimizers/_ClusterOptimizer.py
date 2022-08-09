# -*- coding: utf-8 -*-
"""
Tweaks the cluster centers and spreads, and uses the
CachedClusterAnalyzer class function (Cluster Statistics) to optimize the cluster's 
expected return

Must be passed an instance of CachedClusterAnalyzer, which must be set up
with a datasource


See if for these __init__ statments that if we can assign certain params if specified,
but don't assign if not specified


@author: twbur
@author: srburnett
@author rdhurlbu


Looking at SciPy:
    https://www.w3schools.com/python/scipy_spatial_data.asp
    
    SciPy spacial data seems to have convenient functions that we could use!
        Cosine distance
        Distance
        KDTrees: nearest neighbors!
    
    I think our methods are better than these, as in measuring distance in terms
    of local spreads... It is kinda cool to look at this stuff though
    
    SciPy machine learning
        https://www.w3schools.com/python/python_ml_multiple_regression.asp
        https://www.w3schools.com/python/python_ml_polynomial_regression.asp
        https://www.w3schools.com/python/python_ml_train_test.asp
        
    In general it may be kinda bad to look too much at existing machine learning stuff
    since it may cause you to be less creative
    
    We might actually want to use train and test, that is kinda what we do naturally
    by adjusting stuff and seeing what it does
    
    We definately need to cluster first, since the datapoints are so random.
    However, once we have the clusters, if the clusters follow a pattern like
    reflection over 0, that would suggest that the volitility matters rather
    than the trend, for example. Also, if the clusters draw out a line or a curve
    it would suggest that a combination of dimensions such as addition or 
    subtraction is what matters, rather than the raw dimensions individually
    
    1. Get the clusters
    2. Analyze the total group of clusters for patterns (see if we can refactor a dimension
                                                         to combine them)
    3. Analyze the individual clusters (see if we can add another dimension or change
                                        a dimension to separate the good from the bad)
    
    Definately maybe check out file handling, it would be nice to have ouptut files,
    especially for the stuff that takes forever
    
    Finding clusters is like a shortcut for nearest neighbors. Technically, you could
    go and find the nearest neighbors for each stock on a day, but that would take forever.
    It's easier to find the good areas once, and then just buy a stock if it lands in the
    good area.
    However, when it comes to grouping multiple small clusters together, nearest neighbors
    becomes pretty attractive. We need the clusters first though
        
    Make the cached cluster analyzer be the child of the cluster optimizer?
    They both want to know the cen_spr and the adjust_cen_spr (maybe adj_ratio too)
    
    Have the adj_ratio be a property of the cluster optimizer rather than the cluster analyzer

"""

from copy import deepcopy

#import numpy as np

from _GeneralFunctions import Round_SF

class ClusterOptimizer :
    def __init__(self, clusterAnalyzer, cen_spr, adjust_cen_spr, adj_ratio) :
        self.clusterAnalyzer = clusterAnalyzer
        
        self.adj_ratio = adj_ratio
        
        clusterAnalyzer.cen_spr        = cen_spr
        clusterAnalyzer.adjust_cen_spr = adjust_cen_spr
        
        # Get vars from the clusterAnalyzer
        # TODO make the clusterAnalyzer self tune its adj_ratio or wiggle room
        #adj_ratio = clusterAnalyzer.adj_ratio
        #dataSource = clusterAnalyzer.dataSource
        
        # Spread ratios (for 2: same, halved, doubled)
        # self.spread_adjustments = [1, 1/adj_ratio, adj_ratio]
        #self.spread_adjustments = [1, 1/adj_ratio/adj_ratio, 1/adj_ratio, adj_ratio, adj_ratio * adj_ratio]
        self.spread_adjustments = [1, adj_ratio * adj_ratio, adj_ratio, 1/adj_ratio, 1/adj_ratio/adj_ratio]
        
        # How much to shift a curve center
        delta1 = (adj_ratio - 1)# * dataSource.sqrt_n_dims
        
        delta2 = (adj_ratio * adj_ratio - 1)# * dataSource.sqrt_n_dims
        
        
        # Box mode 
        #if clusterAnalyzer.cache_drop < 1 :
        #    delta *= clusterAnalyzer.cache_drop
        
        # self.center_adjustments = [0, -delta, delta]
        self.center_adjustments = [0, -delta2, delta2, -delta1, delta1]
        
        #                 same        bigger        right         left            smaller
        self.adjustments = [[0, 1], [0, adj_ratio], [delta1, 1], [-delta1, 1], [0, 1/adj_ratio]]
        
        # Save whether the center or spread is adjusted, by index
        self.adjstd_cs  = [[True,  False,          True,        True,         False],  # Center
                           [True,  True,           False,       False,        True]]   # Spread
        
        #self.delta = delta
        
        self.record = []
        
    
    # c_s is an array of the center and spread indexes we will optimize
    # [0, 1] means optimize center first, then spread
    # [1] means optimize spread only
    def Optimize_Step(self, c_s) :
        """
        Try to maximize the Score_Curve function by iterating through all dimensions
        and incrementing in the best direction
    
        Returns
        -------
        None.
    
        """
        
        clusterAnalyzer = self.clusterAnalyzer
        cen_spr = clusterAnalyzer.cen_spr
        adjust_cen_spr = clusterAnalyzer.adjust_cen_spr
        dataSource = clusterAnalyzer.dataSource
        n_dims = dataSource.n_dims
        
        #spread_adjustments = self.spread_adjustments
        #center_adjustments = self.center_adjustments
        adjustments = self.adjustments
        adjstd_cs = self.adjstd_cs
        
        og_set = deepcopy(cen_spr)
        max_score = [-999]
        changed = False
        
        # Hold a record of everything we've tried
        step_record = ["STEP START", og_set]
        
        # Cycle through the dimensions we're gonna optimize
        for dim_idx in range(n_dims):
            
            # Create a record for the current dimension
            step_dim_record = [dataSource.dim_names[dim_idx]]
            
            #keep track of the highest score we get accross this dimension's adjustments
            score = -9999
            max_score = [score]
            max_dim_val = [cen_spr[dim_idx][0], cen_spr[dim_idx][1]]
            
            # Cycle through the adjustments we're gonna make
            for adj_i in range(len(adjustments)) :
                # See if we're gonna make this adjustment
                # Are we adjusting something that impacts this current adjustment?
                if (adjust_cen_spr[dim_idx][0] & adjstd_cs[0][adj_i]) | (adjust_cen_spr[dim_idx][1] & adjstd_cs[1][adj_i]) :
                    # Set the center to the adjusted center
                    cen_spr[dim_idx][0] = og_set[dim_idx][0] + og_set[dim_idx][1] * adjustments[adj_i][0]
                    # Set the spread to the adjusted spread
                    cen_spr[dim_idx][1] = og_set[dim_idx][1] * adjustments[adj_i][1]
                    
                    # If we made the adjustment, get the score
                    score = clusterAnalyzer.Score_Curve(cen_spr)
                    
                    # Record what we got for troubleshooting and future reference
                    step_dim_record.append("[" + (Round_SF(clusterAnalyzer.boundaries[0][dim_idx][0], 4) + ", " +
                                            Round_SF(clusterAnalyzer.boundaries[0][dim_idx][1], 4) + "], " +
                                            Round_SF(score[0], 4)))
                    
                    # Did we get a better score? 
                    if(score[0] > max_score[0]):
                        max_score = score
                        max_dim_val[0] = cen_spr[dim_idx][0]
                        max_dim_val[1] = cen_spr[dim_idx][1]
                        # Assumes that adj = 0 means no change
                        if adj_i:
                            changed = True
            
            # Set the new optimal value to the cen_spr we're tweaking
            cen_spr[dim_idx][0] = max_dim_val[0]
            cen_spr[dim_idx][1] = max_dim_val[1]
                    
            
            """
            # For center and spread 
            for cs in c_s :
                # Is this a center or spread we can adjust?
                if(adjust_cen_spr[dim_idx][cs]):
                    #keep track of the highest total stat score across the 3 adjustments
                    score = -9999
                    max_score = [score]
                    max_dim_val = cen_spr[dim_idx][cs]
                    
                    # If we're optimizing spread
                    if cs :
                        adjustment_range = len(spread_adjustments)
                    # If we're optimizing center
                    else :
                        adjustment_range = len(center_adjustments)
                        
                    # Same, Decrement, Increment
                    for adjustment in range(adjustment_range) :
                        if cs : # Adjust Spread 
                            cen_spr[dim_idx][1] = og_set[dim_idx][1] * spread_adjustments[adjustment]
                        else: # Adjust Center
                            cen_spr[dim_idx][0] = og_set[dim_idx][0] + og_set[dim_idx][1] * center_adjustments[adjustment]
                        
                        #print(cen_spr)
                        #print(cen_spr[dim_idx][cs])
                        
                        score = clusterAnalyzer.Score_Curve(cen_spr)
                        
                        #step_record.append((cen_spr, score[0]))
                        #step_record.append(dataSource.dim_names[dim_idx] + ": " + str(clusterAnalyzer.boundaries[0][dim_idx]))
                        step_dim_record.append("[" + (Round_SF(clusterAnalyzer.boundaries[0][dim_idx][0], 4) + ", " +
                                               Round_SF(clusterAnalyzer.boundaries[0][dim_idx][1], 4) + "], " +
                                               Round_SF(score[0], 4)))
                        #print(score)
                        # Did we get a better score? 
                        if(score[0] > max_score[0]):
                            max_score = score
                            max_dim_val = cen_spr[dim_idx][cs]
                            # Assumes that adj = 0 means no change
                            if adjustment:
                                changed = True
                    # Set the new optimal value to the cen_spr we're tweaking
                    cen_spr[dim_idx][cs] = max_dim_val
                    
                    
                    
                    # Remember the new optimal value
                    #og_set[dim_idx][cs] = max_dim_val
                    
                    # Append the new optimal value to the record
                    #step_record.append((np.array(clusterAnalyzer.boundaries[0]), max_score[0], np.array(clusterAnalyzer.boundaries[1]) - np.array(clusterAnalyzer.boundaries[0])))
                    #step_record.append((cen_spr, max_score))
                    
                    # Stuff with the same cen_spr can get a different (lower) result from round to round,
                    # given that the after-opt boundaries are all included within another
                    # The after-opt boundaries can and should be ok being different from each other from
                    # end opti to another, since different ranges were tried in each end opti
                #print(cen_spr)
            """
            step_record.append(step_dim_record)
        self.record.append(step_record)
        return (changed, max_score)
    
    
    def Optimize_Curve(self, cen_spr) :
        self.cen_spr = cen_spr
        
        print("optimize stat step record")
        stats = 0
        
        # Goal for adj_cen_spr
        # [[],     #Don't optimize center or spread for dim0 (time)
        # [0,1],   #Optimize center and spread for      dim1
        # [1]]     #Optimize spread only for            dim2
        
        while True:
            stats = self.Optimize_Step([0,1])
            #stats = self.Optimize_Step([1])
            #print(stats)
            print()
            self.clusterAnalyzer.dataSource.printCenSpr(cen_spr)
            #print(stats[1])
            print()
            self.clusterAnalyzer.Print_Stats(stats[1])
            print()
            if(not stats[0]) :
                break
    
