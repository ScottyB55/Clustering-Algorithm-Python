# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 00:57:26 2021

@author: srbur
"""

from copy import deepcopy
import math
import matplotlib.pyplot as plt

class BackTester :
    def __init__(self, ClusterOptimizer, Strategy, start_day = 0, skip_initial_opt = True, reopt_period = 9999, start_money = 10000) :
        self.ClusterOptimizer = ClusterOptimizer
        self.ClusterAnalyzer = ClusterOptimizer.clusterAnalyzer
        ClusterAnalyzer = self.ClusterAnalyzer
        
        self.Strategy = Strategy
        
        #self.DataPointSource = ClusterAnalyzer.dataSource
        #DataPointSource = self.DataPointSource
        #market_cap_index = DataPointSource.MarketCap
        
        dataSource = ClusterAnalyzer.dataSource
        cen_spr = ClusterAnalyzer.cen_spr
        #market_cap_index = dataSource.MarketCap
        #ovn_ret_index = dataSource.OvernightRet
        
        # Save the number of days that we have gone without reoptimizing
        reopt_count = 0
        
        # How many days we go before reoptimizing
        #reopt_period = 99999#62
        # Is the initial thing optimized (Should we skip optimization at beginning?)
        #optimized = False#True
        optimized = skip_initial_opt
        
        self.picks_by_day = []
        
        self.start_day = start_day
        
        current_day_on = start_day
        
        
        days_points = []
        #self.days_picks = days_picks
        
        
        #days_count = 0
        
        #days_weight = 0
        
        #days_returns = []
        
        #returns_by_day = []
        
        #total_return = 0
        
        total_count = 0
        
        
        
        total_money = start_money
        # Start off with 10k
        #total_money = 10000
        #total_money = 1000
        #total_money = 1
        #total_money = 0.001
        #total_money = 1000000
        
        # We subtract 10 from the logged market cap (divide by about 20k, 0.005%)
        #market_fraction = -10
        
        # Keep track of the total moneys across time
        # Don't buy more than the specified fraction of the market cap for any stock
        self.money_by_day = [total_money]#[]
        
        # Datapoint structure:
        # Day is index 0, market cap is index 1, overnight return (dep var) is index n_dims + 1 (n_dims is watered down low risk)
        
        # Make sure the cache is up to date around the current cen_spr
        ClusterAnalyzer.Update_Parent_Cache(start_day, 9999, True, True)
        ClusterAnalyzer.Update_Cache(cen_spr)
        
        # Add on a fake point with the next day to cause the loop to loop one more time
        fake_point = deepcopy(ClusterAnalyzer.caches[0][-1])
        fake_point[0] += 1
        
        for point in (ClusterAnalyzer.caches[0] + [fake_point]) :
            # Only do the days after we started
            if point[0] >= start_day :
                if (not optimized) & (not (ClusterOptimizer is None)) :
                    # Allow the cluster analyzer to see data up until the current day (exclusive)
                    cen_spr[0][0] = current_day_on
                    #ClusterAnalyzer.Update_Time_Range(0, current_day_on, True)
                    ClusterAnalyzer.Update_Parent_Cache(0, current_day_on, True, True)
                    # Optimize the cluster
                    # TODO make sure that the center spread that the optimizer starts with 
                    # is the most recently updated one! Make sure that if the centers and 
                    # spreads are saved in many different places, that they all point to
                    # the same location!
                    ClusterOptimizer.Optimize_Curve(cen_spr)
                    optimized = True
                    
                    print("Optimized up to day " + str(current_day_on) + "!")
                    print("Cluster center and spread:")
                    #print(cen_spr)
                    
                    dataSource.printCenSprArray(cen_spr)
                    print()
                    #dataSource.printCenSpr(cen_spr)
                    #print()
                    
                # If we moved up a day, add the previous day's array and reset for the new day
                while point[0] > current_day_on :
                    
                    # repeatedly cycle through the days_picks that we haven't bought yet
                    # Save the initial stocks, and money we have to buy
                    
                    
                    days_buys = self.Strategy.Get_Buys(days_points, total_money)
                    
                    days_picks = days_buys[0]
                    leftover_money = days_buys[1]
                    
                    days_buys = self.Strategy.Get_Sales(days_picks, leftover_money)
                    
                    days_picks = days_buys[0]
                    total_money = days_buys[1]
                    
                    self.money_by_day.append(total_money)
                    
                    # Save the day's picks into an array of the picks for each day
                    self.picks_by_day.append(days_picks)#deepcopy(days_picks)) we redefine the pointer below, deepcopy not necessary
                    
                    total_count += len(days_picks)
                    
                    # Cycle through the stuff we bought today
                    # Find the average return over the day
                    
                    
                    self.days_picks = days_picks
                    days_points = []
                    
                    
                    #days_count = 0
                    
                    #days_weight = 0
                    
                    current_day_on += 1
                    
                    # Increment the day count until we reoptimize
                    reopt_count += 1
                    
                    # If we hit the day count for reoptimization, enable optimizaion
                    if reopt_count >= reopt_period :
                        reopt_count = 0
                        optimized = False
                
                # See if the current point fits the curve
                #weight = ClusterAnalyzer.Get_Weight(point)
                
                days_points.append(point)
                
                    
                
        
        print("Total Count")
        print(total_count)
        
        #We need to compare the squared distance to cache_drop squared rather than cache drop!
        
        #Plot it!
        
        # TODO make the plot x axis reflect the range from the start date to current day on
        # ValueError: x and y must have same first dimension, but have shapes (304,) and (184,)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        #ax.set_yscale('log')
        #ax.plot(range(current_day_on-start_day), returns_by_day)
        #ax.plot(range(current_day_on-start_day), money_by_day)
        ax.plot(range(current_day_on-start_day+1), self.money_by_day)
        plt.xticks(rotation=45)
        plt.title('Growth of 10k, factoring in saturation')
        
        # We need to factor in the number of days that we have points when rating the curve
        # Average return * confidence * (# of days active / total # of days)
        # Get a weighted number of days somehow