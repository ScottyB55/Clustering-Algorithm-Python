# -*- coding: utf-8 -*-
"""
Takes a data source and a cluster we're analyzing (set of centers and spreads)

Analyzes the cluster, and fits points to the cluster (filters out insignificant points from its cache memory)

Cache eliminates excessive looping through insignificant data points

Wiggle room in the cache can allow for adjustments without having to update the cache

The cache currently appears to have a slight problem about sometimes not containing all the significant datapoints

**Let's set the cache standard to the trapezoid cluster shape!!


Function outputs:
    Datapoint cluster fit, Point Weight
    Cluster Score, Cluster Statistics

When the cluster or datasource is updated:
    When anything is updated (cluster or datasource), we must validate or update the cache
    (ensure that the cache contains all relevant datapoints!)


Smart cache is adaptive and can be used by any optimization algorithm
It can be used by any optimization algorithm


@author: twbur
@author: srburnett
@author rdhurlbu
"""

import math
from copy import deepcopy

from _GeneralFunctions import Round_SF

class CachedClusterAnalyzer :
    
    # Initialization function!
    def __init__(self, dataSource, cen_spr, day_range = (0,9999), strat_tradeoff = 180, point_n_confidence = 40, cluster_type = None, ln_daily_n = True) :        
        
        
        # What day range of dataSource.dataPoints are we looking at
        self.min_day = day_range[0]   #day_range = (min_day, max_day)
        self.max_day = day_range[1]
        
        if cluster_type == "sloped box" :
            self.box_mode = True
            self.gradual_dropoff = True
            self.quartic_dist = True
        elif cluster_type == "quartic":
            self.box_mode = False
            self.gradual_dropoff = True
            self.quartic_dist = True
        else :
            if not (cluster_type == "steep curve") :
                print("Cluster type assigned to steep curve by default")
            self.box_mode = False
            self.gradual_dropoff = False
            self.quartic_dist = False
        
        
        self.dataSource = dataSource
        
        self.cen_spr = cen_spr
        
        # This allows us to scale down the weighted number of points in a day
        # We can have two options: ln and square root
        self.scale_down = []
        for i in range(8000) :
            if i :
                if ln_daily_n :
                    self.scale_down.append(math.log(i+1) / i)
                else :
                    self.scale_down.append(math.sqrt(i) / i)
            else :
                self.scale_down.append(1)
        
        
        
        # We will hold the starting index of each day in dataSource.dataPoints, easy day start lookup
        self.day_keys = []
        day_keys = self.day_keys
        
        # Fill out the day_keys array!
        current_day=0
        current_index=0
        current_day_start_index=0
        
        # Cycle through all the dataPoints in the dataSource (Must be sorted by day!)
        # TODO we already do this in the new datasource, so let's eventually use the
        # one in there instead of being redundant!
        for point in dataSource.DataPoints :
            # While the incoming day is greater than the current day,
            # append the current day start index, increment the day, 
            # and set the next current day start index to the next index
            while point[0] > current_day :
                day_keys.append(current_day_start_index)
                current_day += 1
                current_day_start_index = current_index# + 1
            
            # increment the current index
            current_index += 1
            
        # add on the last day start index
        day_keys.append(current_day_start_index)
        
        
        # Tradeoff between number of points, variation, and average return  
        self.strat_tradeoff = strat_tradeoff
        
        self.point_n_confidence = point_n_confidence
        
        #%%
        
        
        # Total exponential distance input after which datapoints are insignificant, doesn't matter for box
        self.cache_dropoff = 4
        
        # Box mode generally should use greater adjustment values
        if self.box_mode :
            # Total distance until cutoff, in terms of spreads
            # Means that a distance of 1 spread yields a weight of 0.5 for the trapezoid
            self.dist_2_cutoff = 1.333#dataSource.sqrt_n_dims
            
            # 1 dimensional distance until cutoff (for cache filtering), in terms of spreads
            self.cache_1D_dropoff = self.dist_2_cutoff
            
        else :
        # Curve mode, one spread in each direction yields a weight of e^-1*dataSource.sqrt_n_dims
            if self.gradual_dropoff :
                self.dist_2_cutoff = self.cache_dropoff
            else:
                self.dist_2_cutoff = math.sqrt(self.cache_dropoff)
            self.cache_1D_dropoff = math.sqrt(self.dist_2_cutoff * dataSource.sqrt_n_dims) # * dataSource.n_dims
        
            if self.quartic_dist :
                self.cache_1D_dropoff = math.sqrt(self.cache_1D_dropoff)
        
        # Filter down the datapoints into smaller relevant sections
        #     small    medium    large(parent)
        self.caches = [[], [], []]
        
        # Setup the parent cache to the time filtered datapoints
        # Filter Nans
        self.Update_Parent_Cache(self.min_day, self.max_day, True, True)
        
        
        """
        Boundaries:
        Points of each dimension close to the current center spread, smallest cache, and medium cache
        Outer boundary (what range does the group contain)
        """
        
        # How many updates to a boundary do we save? (the history allows us to calculate the parent boundary)
        #self.num_hist  = [5, 25, 1]#* 7   # 3 * 2 * 6
        self.num_hist  = [3, 24, 1]
        # Save the history index we're on for each boundary
        self.bound_hist_ind = [0] * len(self.caches)
        
        # Define the current range of values we're looking at across each dimension, for each boundary
        self.boundaries = [[ [0.0] * 2 for _d in range(dataSource.n_dims)] for _c in range(len(self.caches))]
        
        # Save the history of the boundaries array, that we analyze to calculate new boundaries
        self.boundary_history =     [[[ [999999.9, -999999.9] for _d in range(dataSource.n_dims)] for _t in range(self.num_hist[_c])] for _c in range(len(self.caches))]
        
        # Define a net range measurement that we will make on the boundary_history
        self.history_range =       [[ [0.0] * 2 for _d in range(dataSource.n_dims)] for _c in range(len(self.caches))]
        
        #%%



    #%%
    def Update_Cache(self, cen_spr):
        # Get the self values
        dataSource = self.dataSource
        n_dims = dataSource.n_dims
        self.cen_spr = cen_spr
        #cache_drop_1D = self.cache_drop_1D
        cache_1D_dropoff = self.cache_1D_dropoff
        #outer_boundary = self.outer_boundary
        #inner_boundary = self.inner_boundary
        #cache_buffer = self.cache_buffer
        
        caches = self.caches
        
        boundaries_Updated = self.boundaries_Updated
        
        num_hist = self.num_hist
        bound_hist_ind = self.bound_hist_ind
        boundaries = self.boundaries
        boundary_history = self.boundary_history
        history_range = self.history_range
        
        """ Called when the center or spread is changed. Will update the cache when 
        the cache needs updating
        
        Finds the caches that need to be updated then updates them. 
        
        P
        
        """
        
    #%%    
        # Update the current boundary based on the center and spread
        for i in range(n_dims) :            
            center = cen_spr[i][0]
            spread = cen_spr[i][1] * cache_1D_dropoff
            
            boundaries[0][i][0] = center - spread
            boundaries[0][i][1] = center + spread
        
        # Reset the old_bound index if we spill over
        bound_hist_ind[0] += 1
        if(bound_hist_ind[0] >= num_hist[0]) :
            bound_hist_ind[0] = 0
        
        
        # Which boundaries were updated?
        # We updated the smallest boundary
        #boundaries_Updated = [True, False, False]
        boundaries_Updated[0] = True
        
        # Cycle through boundaries 
        for boundary in range(0, len(caches)-1):
            # Update the boundary history with the current boundary
            boundary_history[boundary][bound_hist_ind[boundary]] = deepcopy(boundaries[boundary])
            
            # Keep tabs on the boundary's history range
            bound_hist_range       = history_range[boundary]
            
            
            # Keep a tab on the current boundary
            current_boundary = boundaries[boundary]#boundary_history[boundary][bound_hist_ind[boundary]]
            
            # Keep tabs on the container boundary (to see if the current boundary crosses it)
            container_boundary = boundaries[boundary+1]
            
            
            # If the boundary was updated, calculate its statistics!
            if(boundaries_Updated[boundary]) :
                
                # Calculate statistics on the boundary's history
                # TODO we only need to do this if the current boundary changes!
                # (we are finding the range, or min and max for now, for each dimension)
                for i in range(n_dims) :
                    # Find the historic min and max of this dimension
                    net_min = 999999.9
                    net_max = -999999.9
                    for old_bound in boundary_history[boundary]:
                        if old_bound[i][0] < net_min :
                            net_min = old_bound[i][0]
                            #print(net_min)
                        if old_bound[i][1] > net_max :
                            net_max = old_bound[i][1]
                            
                    # Update the boundary's historic range
                    bound_hist_range[i][0] = net_min
                    bound_hist_range[i][1] = net_max
                
                for i in range(n_dims) :
                    # See if the current (enclosed) boundary exceeds the boundary that contains it
                    if          (current_boundary[i][0] < container_boundary[i][0]
                            ) | (current_boundary[i][1] > container_boundary[i][1]) :
                        # If the enclosed boundary exceeds the boundary that contains it,
                        # We need to update the container boundary!
                        boundaries_Updated[boundary+1] = True
                        break
                
                
                
                # If the outer boundary (container_boundary) needs updating, update as a fn of:
                    # Current boundary (current_boundary)
                    # Current net statistics (bound_hist_range)
                if boundaries_Updated[boundary+1] :
                    # TODO potentially make this boundary update function more complex
                    
                    boundaries[boundary+1] = deepcopy(bound_hist_range)
                    
                    # Update the boundary history index, and reset if it spills over!
                    # Reset the old_bound index if we spill over
                    bound_hist_ind[boundary+1] += 1
                    if(bound_hist_ind[boundary+1] >= num_hist[boundary+1]) :
                        bound_hist_ind[boundary+1] = 0
                
                
                
                
                #%%
        """ Update the caches now that we know their boundaries 
        The small caches draw from the big caches, so we update the big bois first"""
        # For each cache, starting from the largest mutable one
        for boundary in range(len(caches) - 1, 0, -1):
            
            #see if this cache needs updating
            if (boundaries_Updated[boundary]):
                # We updated the cache to the new boundary, so set it to non-updated for the future
                boundaries_Updated[boundary] = False
                # Reset the cache
                caches[boundary-1] = []
                # cycle through the parent cache
                #print(boundaries[boundary][0])
                for point in caches[boundary] :
                    fits = True
                    for i in range(n_dims) :
                        if i :
                            dimension = point[i]
                            # The boundary index is the corresponding cache index + 1!
                            dimension_min = boundaries[boundary][i][0]
                            dimension_max = boundaries[boundary][i][1]
        
                            if (dimension < dimension_min) | (dimension > dimension_max):
                                fits = False
                                break
                    #if (point[0] < 308) & (point[0] > 280) :
                    #    print (point[0])
                    if fits :
                        #print(point[0])
                        caches[boundary-1].append(point)
                print("Updating cache " + str(boundary-1) + ", " + str(bound_hist_ind[boundary-1]) + ", " + str(len(caches[boundary-1])))
    #%%
    
    
    
    
    
    
    
    
    """ THIS IS THE END OF THE CACHE SECTION
    
        BEGINNING THE CURVE ANALYSIS SECTION
    """
    
    
    
    
    
    
    
    
    #%%
    # Calculate the weight given cen_spr for a particular point in the point_cache, indexed by p_c_index
    # For boxes and curves, depending on cache_drop!
    def Get_Weight(self, point):
        # Get the self values
        dataSource = self.dataSource
        n_dims = dataSource.n_dims
        sqrt_n_dims = dataSource.sqrt_n_dims
        quartic_dist = self.quartic_dist
        cen_spr = self.cen_spr
        #cache_drop = self.cache_drop
        #cache_drop2 = self.cache_drop2
        dist_2_cutoff = self.dist_2_cutoff
        gradual_dropoff = self.gradual_dropoff
        box_mode = self.box_mode        
        
        
        """ Gets a weight on the given curve for a single point
    
        Parameters
        ----------
        point : The point to weigh
    
        Returns
        -------
        If we are in box mode, we will only return a value of 1 or 0, if they are within the boundary
        If we are in curve mode, we will return a continueous value for weight
        
        """
        
        # Box mode!
        if box_mode :
            # Do the hard cutoff if we have the extra dropoff
            if not gradual_dropoff :
                for i in range(n_dims) :
                    # get the center and spread for the current dimension
                    cen = cen_spr[i][0]
                    spr = cen_spr[i][1] * dist_2_cutoff
                    
                    # get the point's location along the current dimension
                    xvar = point[i]
                    
                    # calculate the distance between the point and center
                    dist = xvar - cen
                    
                    # if we are out of range, return 0 weight
                    if((dist > spr)) | ((dist < -spr)) :
                        return 0
                return 1
            # Do the ramped cutoff if we want a slower dropoff
            # save the lowest ramp or value we hit
            ramp_min = 1
            for i in range(n_dims) :
                # get the center and spread for the current dimension
                cen = cen_spr[i][0]
                spr = cen_spr[i][1] * dist_2_cutoff
                
                # get the point's location along the current dimension
                xvar = point[i]
                
                # calculate the distance between the point and center
                dist = xvar - cen
                
                # get the absolute value of the distance
                if (dist < 0) :
                    dist = -dist
                
                # if we are out of range, return 0 weight
                if (dist >= spr) :
                    return 0
                
                # if we are completley in range, don't change anything
                dist += dist # multiply distance by 2
                
                if (dist <= spr) :
                    continue
                    
                # if we are on the ramp, update the ramp_min
                dist = 2 - (dist / spr)
                if (dist < ramp_min) :
                    ramp_min = dist
                
            return ramp_min
        
        # Curve Mode!
        cen_dist = 0
        for i in range(n_dims) :
            # get the center and spread for the current dimension
            cen = cen_spr[i][0]
            spr = cen_spr[i][1]
            
            # get the point's location along the current dimension
            xvar = point[i]
            
            # calculate the distance between the point and center in terms of the spread
            dist = (xvar - cen) / spr
            
            if quartic_dist :
                dist *= dist
            
            # factor this into the total distance
            cen_dist += dist * dist / sqrt_n_dims#n_dims
            
            # if we are out of range, return 0 weight
            if (cen_dist > dist_2_cutoff):
                return 0
        if (gradual_dropoff) :
            return math.exp(-cen_dist)
        return math.exp(-cen_dist*cen_dist)
        
    #%%
    
    def Score_Curve(self, cen_spr) :
        # Get the self values
        self.cen_spr = cen_spr
        dataSource = self.dataSource
        n_dims = dataSource.n_dims
        caches = self.caches
        strat_tradeoff = self.strat_tradeoff
        point_n_confidence = self.point_n_confidence
        
        """
        Calculates statistics about a given curve.
    
        Returns
        -------
        A tuple representing aspects of our curve. This is in the form: 
             (curve_score, average_ret, average_var, total_weighted_n, effective_n, confidence)
            
        """
        
        self.Update_Cache(cen_spr)
        total_weighted_n = 0
        total_weighted_ret = 0
        total_weighted_var = 0
        
        day_weighted_n = 0
        day_weighted_ret = 0
        day_weighted_var = 0
        
        # The cache is organized by day!
        # Let's cycle through the cache and find the number of weighted points per day
        
        # Save the day we're on
        current_day_on = 0
        
        # Which cache we're gonna use. Used in testing to see if there's a difference between the caches.
        # They should all yield the same result, the only difference being that using the larger cache
        # Should be slower
        
        """
        cache_num = 1
        
        Cluster Score: 0.5953
        Average Return: 0.07236
        Average Variation: 0.08598
        Total Weighted N: 122.9
        
        There's a difference.
        
        
        cache_num = 0
        
        Cluster Score: 0.5855
        Average Return: 0.06949
        Average Variation: 0.08187
        Total Weighted N: 122.8
        
        Cluster Score: 0.535
        Average Return: 0.05912
        Average Variation: 0.07259
        Total Weighted N: 144.8
        
        
        We found at least part of the problem! We weren't counting the points where we moved onto a new day!
        
        cache_num = 1
        
        Cluster Score: 0.6119
        Average Return: 0.07326
        Average Variation: 0.08681
        Total Weighted N: 123.7
        
        ITS DOING THE SAME THING!! IT DONE FIXED!! YAAY!!
                
        Cluster Score: 0.6119
        Average Return: 0.07326
        Average Variation: 0.08681
        Total Weighted N: 123.7
        
        SWITCHING BETWEEN THE CACHES USED TO CALCULATE THE SCORE IS A GREAT WAY TO TEST THE CACHE!
        
        cache_num = 0 and cache_num = 1 get the same thing here!
        Cluster Score: 0.3169
        Average Return: 0.05064
        Average Variation: 0.06197
        Total Weighted N: 166.1
        
        
        """
        cache_num = 0
        
        # Cycle through the cache
        for i in range(len(caches[cache_num])) :
            # Get the point of the cache we're on
            point = caches[cache_num][i]
            
            # if we moved onto a new day
            if not (point[0] == current_day_on) :
                # Do operations on the old day's stats
                # Basically scale it down if there's a ton in one day
                #scale_down = math.log(total_weighted_n + 1)
                
                # The scale down is really slowing it down!
                # We need to create a search index plus slope!
                
                if day_weighted_n :
                    scale_d = self.scale_down[round(day_weighted_n)]
                    
                    day_weighted_n *= scale_d
                    day_weighted_ret *= scale_d
                    day_weighted_var *= scale_d
                
                    # Push the old day's stats into the total
                    total_weighted_n += day_weighted_n
                    total_weighted_ret += day_weighted_ret
                    total_weighted_var += day_weighted_var
                    
                    # Reset the day's stats
                    day_weighted_n = 0
                    day_weighted_ret = 0
                    day_weighted_var = 0
                
                # Update the day we're on
                current_day_on = point[0]
            
            # Add the point's weight, weighted return, and variation to the day's culmulative
            point_weight = self.Get_Weight(point)
            if point_weight :
                # We're gonna use the low risk return here!
                point_returns = point[n_dims]
                point_weighted_ret = point_weight * point_returns
                day_weighted_n   += point_weight
                day_weighted_ret += point_weighted_ret
                day_weighted_var += abs(point_weighted_ret)
            
            
        
        # TODO this leaves off the most recent day, above
        
        # Sample size is too small
        if (total_weighted_n <= 1) :
            return (-1, -1, -1, total_weighted_n, -1, -1)
        
        # Used to calculate effective n, in case ret is negative
        #abs_total_weighted_ret = abs(total_weighted_ret)
        
        #effective_n = abs_total_weighted_ret / (total_weighted_var - 0.5 * abs_total_weighted_ret) * total_weighted_n
        
        # Tradeoff between datapoints and average return, and some variation            
        # e^-0 = 1 = complete confidence
        #strat_tradeoff = -200
        
        # Approaches e^-0 -> 1 as the effective number of points increases
        #confidence = math.exp( strat_tradeoff / effective_n) * effective_n * self.scale_down[round(effective_n)]
        #confidence = math.exp( strat_tradeoff / effective_n)
        #confidence = math.exp( strat_tradeoff / effective_n) * math.sqrt(effective_n)
        
        # Calculate the average return in the cluster
        average_ret = total_weighted_ret / total_weighted_n
        abs_average_ret = abs(average_ret)
        
        # Calculate the average variation in the cluster
        average_var = total_weighted_var / total_weighted_n
        
        # Calculate the confidence ratio (that its not random)
        # Actually, now we're using this as a way to try to get a datapoint a day
        ret_var_ratio = abs_average_ret / (average_var - 0.5 * abs_average_ret)
        
        effective_n = total_weighted_n * ret_var_ratio
        
        tradeoff_input = strat_tradeoff / (effective_n)
        tradeoff_input *= tradeoff_input
        
        # Factors in the confidence that it's not random
        conf_input = point_n_confidence / (effective_n)
        
        confidence = math.exp ( - conf_input - tradeoff_input )
        
        # Calculate the confidence that it's not random
        #conf_input = strat_tradeoff / (total_weighted_n * ret_var_ratio)
        #confidence = math.exp( conf_input )
        
        
        #curve_score = average_ret * confidence
        curve_score = average_ret * ret_var_ratio * confidence# * math.sqrt(total_weighted_n)
       
        #return (curve_score, average_ret, average_var, total_weighted_n, effective_n, confidence)
        return (curve_score, average_ret, average_var, total_weighted_n, ret_var_ratio, confidence)
    
    def Print_Stats(self, stats) :
        print("Cluster Score: " + Round_SF(stats[0], 4))
        print("Average Return: " + Round_SF(stats[1], 4))
        print("Average Variation: " + Round_SF(stats[2], 4))
        print("Total Weighted N: " + Round_SF(stats[3], 4))
    
    def Update_Parent_Cache(self, min_day, max_day, bool_min_update, bool_nan_filter) :
        self.min_day = min_day
        self.max_day = max_day
        
        # Update the min time cutoff if we enabled it
        if bool_min_update :
            self.Update_Min_Time_Cutoff()
            min_day = self.min_day
            max_day = self.max_day
        # TODO: verify the following in all functions
        # IF WE CALL A FUNCTION THAT UPDATES THE SELF.VALUES,
        # WE MUST DEFINE THE SELF.VALUES AFTER THAT FUNCTION!
            
        print("min day: " + str(min_day))
        print("min index: " + str(self.day_keys[min_day]))
        
        # TODO make an optional nan filter too!
        
        # Filter the datapoints down by the new time range
        # This is actually by value, the colon returns a copy of the range!
        if max_day >= len(self.day_keys) :
            timeFilteredDataPoints = self.dataSource.DataPoints[self.day_keys[min_day] : ]
        else :
            timeFilteredDataPoints = self.dataSource.DataPoints[self.day_keys[min_day] : self.day_keys[max_day]]
        
        # Since the colon returns a copied value, we can delete indexes
        # Without affecting the origional datasource!
        if bool_nan_filter :
            FilteredDataPoints = []
            for datapoint in timeFilteredDataPoints :
                add_point = True
                for var in datapoint :
                    if math.isnan(var) :
                        add_point = False
                        break
                if add_point :
                    FilteredDataPoints.append(datapoint)
            self.caches[len(self.caches)-1] = FilteredDataPoints
        else :
            self.caches[len(self.caches)-1] = timeFilteredDataPoints
                
            """
            for datapoint in self.FilteredDataPoints :
                for var in datapoint :
                    if math.isnan(var) :
                        del datapoint
                        break
            """
        
        # Update the parent cache to the new time range
        # The highest level of cache is a reference to timeFilteredDataPoints
        #self.caches[len(self.caches)-1] = self.FilteredDataPoints
        
        print("Parent cache (filtered) len: " + str(len(self.caches[len(self.caches)-1])))
        
        # Show that we need to update all caches since there is new data!
        self.boundaries_Updated = [True] * len(self.caches)
    
    # Updates the min_day to the time we can cut off at, given the time center, spread, and dropoff
    def Update_Min_Time_Cutoff(self) :
        min_day = self.min_day
        #max_day = self.max_day
        cen_spr = self.cen_spr
        cache_1D_dropoff = self.cache_1D_dropoff
        
        # The first cen_spr is reserved for time,
        # Let's see how many spreads back we need to go to dropoff
        drop_day = cen_spr[0][0] - cen_spr[0][1] * cache_1D_dropoff
        
        drop_day = round(drop_day)
        
        # Make sure its not negative
        if drop_day < 0 :
            drop_day = 0
        
        if min_day < drop_day :
            self.min_day = drop_day