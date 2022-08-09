# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 17:21:12 2021

@author: srbur
"""

from copy import deepcopy
import math

class Strategy :
    
    # TODO make the max fraction be a function of the point's weight! Done:
    # money_limit = min(money_saturate, max_money_per) * math.sqrt(day_point_info["Weight"])
    def __init__(self, ClusterAnalyzer, buy_threshold = 0.1, buy_limit = 10, max_fraction = 0.33, market_fraction = -9, milk_curve = True) :
        self.ClusterAnalyzer = ClusterAnalyzer
        self.buy_threshold = buy_threshold
        self.buy_limit = buy_limit
        self.max_fraction = max_fraction
        # We subtract 10 from the logged market cap (divide by about 20k, 0.005%)
        self.market_fraction = market_fraction
        self.milk_curve = milk_curve
        
        self.dataSource = ClusterAnalyzer.dataSource
        #cen_spr = ClusterAnalyzer.cen_spr
        self.market_cap_index = self.dataSource.MarketCap
        self.ovn_ret_index = self.dataSource.OvernightRet
        
    """ Pass it an array of datapoints (should all be on the same day)
        that we're thinking about buying. And how much money we have
        
        It returns an array of tuples: [(datapoint_to_buy, dollar_value_to_buy),...]
        
        It returns a tuple consisting of:
            a list of dictionaries: days_picks {Index, Point, Money In, Money Out (nan since we don't sell yet), weight}
            the leftover money we (should) have after buying all the stonkies
    """
    def Get_Buys(self, days_points_array, total_money) :
        
        days_picks = []
        
        max_money_per = total_money * self.max_fraction
        
        # cycle through the day's points
        for point in days_points_array :
            # See if the current point fits the curve
            weight = self.ClusterAnalyzer.Get_Weight(point)
        
            # Add the point to the day's array if it fits
            # This isn't the final days picks, we pick the top 10 for example!
            if (weight > self.buy_threshold) :
                point_index = self.dataSource.Get_Point_Index(point)
                point_ticker = self.dataSource.Get_Ticker(point_index)
                
                days_picks.append({"Index"     : 0,    
                                   "Point"     : point, #deepcopy not necessary, not changing vals
                                   "Money In"  : 0.0,  
                                   "Money Out" : math.nan,
                                   "Weight"    : weight,
                                   "Ticker"    : point_ticker})
            
        def sort_weight(point): return point["Weight"]
        days_picks.sort(reverse = True, key = sort_weight)
        top_ten_picks = []
        
        for top_ten_index in range(min(self.buy_limit, len(days_picks))) :
            days_picks[top_ten_index]["Index"] = top_ten_index
            top_ten_picks.append(days_picks[top_ten_index])
        
        days_picks = top_ten_picks
        
        leftover_stocks = deepcopy(days_picks)#deepcopy might be necessary, since we're removing vals
        leftover_money = total_money
        
        
        
        while len(leftover_stocks) :
            # Find out how much money we want to put in per stock
            leftover_weight = 0
            money_per = 0
            if self.milk_curve :
                for day_point_info in leftover_stocks :
                    weight = day_point_info["Weight"]
                    leftover_weight += weight
            else :
                money_per = leftover_money / len(leftover_stocks)
            
            #print("leftover_weight: " + str(leftover_weight))
            
            # Save how much money we'd like to put in each stock (not final money in)
            for day_point_info in leftover_stocks :
                if self.milk_curve :
                    weight = day_point_info["Weight"]
                    money_in = leftover_money * weight / leftover_weight
                else :
                    money_in = money_per
                
                days_picks[day_point_info["Index"]]["Money In"] = money_in
                day_point_info["Money In"] = money_in
                # If we're iterating through an object in a for loop,
                # We actually receive the pointer! Rather than a copy!
                # However, if we're iterating a value such as an int we get the raw
                
                #... actually i think that may have woked because we had the line above... idk
                # We definately would need to change it in the days picks array 
                # because the leftover stocks array is a deepcopy of the original days
                # picks... lol, it has the right weight and index but doesn't know the money in
                # ...Oh! if we want to keep the money in for the leftover stocks,
                # we need that line! Turns out day_point_info actually is a pointer to 
                # the leftover stocks array!
                
            #input("Press Enter to continue...")
            # Will we buy all the stocks at this rate? Or do some saturate?
            buy_all = True
            
            #"""
            # See if the market cap allows us to put that much money in
            for day_point_info in leftover_stocks:
                #day_point = day_point_info[1]
                
                money_saturate = math.exp(day_point_info["Point"][self.market_cap_index] + self.market_fraction)
                
                #max_money_risk = max_money_per * math.sqrt(day_point_info["Weight"])
                
                #money_limit = min(money_saturate, max_money_risk)
                
                money_limit = min(money_saturate, max_money_per) * math.sqrt(day_point_info["Weight"])
                
                # TODO also put in a risk management factor, for example that we won't
                # put more than 33% of our money into one stock (we should make it a function
                # of the weight too)
                
                # Se if we would saturate this stock
                if day_point_info["Money In"] > money_limit:
                    # This stock would over saturate if we bought as much
                    # as we wanted
                    # Buy the saturation limit, then recalculate leftover money
                    leftover_money -= money_limit
                    
                    #print(day_point)
                    # Save how much money we put in
                    days_picks[day_point_info["Index"]]["Money In"] = money_limit
                    
                    buy_all = False
                    
                    leftover_stocks.remove(day_point_info)
                    
                    break
            
            
            # If we wouldn't saturate a stock, buy all of them!
            if buy_all :
                for day_point_info in leftover_stocks:
                    #day_point = day_point_info[0]
                    
                    leftover_money -= day_point_info["Money In"]
                    #leftover_money -= money_per
                    
                    #days_picks[day_point_info["Index"]]["Money In"] = money_per
                    
                    # Technically, we could not remove this from the leftover
                    # Stocks and simply break at the end of this loop
                    # Just doing this to make sure everything works as expected
                    leftover_stocks.remove(day_point_info)
        
        return (days_picks, leftover_money)
    
    
    def Get_Sales(self, days_picks_array, leftover_money) :
        
        # Now that we know which stocks we bought and how many we bought,
        # Save how much money we got out of them!... Actually,
        # We probably want to make this separate from the get buys function
        
        for day_point_info in days_picks_array:
            # Calculate and save how much money we got out of the stock (sold it for)
            money_out = day_point_info["Money In"] * math.exp(day_point_info["Point"][self.ovn_ret_index]) 
            #days_picks[day_point_info["Index"]]["Money Out"] = money_out
            # Byref since we're iterating objects!
            day_point_info["Money Out"] = money_out
            
            leftover_money += money_out
        
        return (days_picks_array, leftover_money)
        