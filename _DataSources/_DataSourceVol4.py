# -*- coding: utf-8 -*-
"""
Input Data Configuration File

Sets up a data source
Used directly by the CachedClusterAnalyzer


@author: twbur
@author: srburnett
@author rdhurlbu

Uses the new datafiles that ryan has a script to automatically update pickles

Memory problem with the new datasource!


import gc

gc.collect()

Run after marking what i want to discard with del

del my_object
"""

file_path = 'data-files/'#../
file_ext = '.pkl'

fund_path = 'data-files\Fundamental_info_union.pkl'

import numpy as np
import pandas
from datetime import datetime, timedelta

#from _GeneralFunctions import Round_SF

import math
#from copy import deepcopy


class FileMama :
    def __init__(self, start_day_index) :
        
        """ All of the self attributes are numpy arrays """
        """ Read all the raw data pickles! """
        
        
        
        
        #file_path = '../data-files/'
        #file_ext = '.pkl'
        
        #Purpose: analyze relationship between RH positioning and equity returns. Use pct. return, price, volume as
        #proxies for popular RH stocks.
        
        # Cleaned, better commented version of Analysis1 code
        
        # Import packages and data
        
        
        #%% Get 4:00 and 5 minute dataframes
        
        data = pandas.read_pickle(file_path + 'rh_volume_data' + file_ext)
        data_5min = pandas.read_pickle(file_path + 'rh_5min_data' + file_ext)
        
        
        #%% Convert the dataframes into a usable form
        data_volume = data_5min.loc[:, data_5min.columns.get_level_values(1) == 'volume']
        data_volume_sum = data_volume.sum(level=[0])
        set_missing = set( [i[0] for i in data.columns] ) - set( [i[0] for i in data_volume_sum.columns] )
        set_union = np.sort(list( set([i[0] for i in data.columns]) - set_missing ))
        
        volume_5min_df = data_volume_sum[ set_union ]
        
        
        #%% Extract the dataframe data into arrays
        
        #Import Return Data
        Open_df = data.loc[:, data.columns.get_level_values(1) == 'open'][ set_union ]
        Close_df_0 = data.loc[:, data.columns.get_level_values(1) == 'market_cap'][ set_union ].values / data.loc[:, data.columns.get_level_values(1) == 'shares_outstanding'][ set_union ].values
        Close_df = pandas.DataFrame(Close_df_0)
        
        #Import Volume Data
        Volume_df = data.loc[:, data.columns.get_level_values(1) == 'volume'][ set_union ]
        #Import Avg Volume Data
        Avg_Volume_df = data.loc[:, data.columns.get_level_values(1) == 'average_volume'][ set_union ]
        #Import High Data
        High_df = data.loc[:, data.columns.get_level_values(1) == 'high'][ set_union ]
        # Import Low Data
        Low_df = data.loc[:, data.columns.get_level_values(1) == 'low'][ set_union ]
        
        # Import Fundamental Data
        mkt_cap_df = data.loc[:, data.columns.get_level_values(1) == 'market_cap'][ set_union ]
        
        
        
        
        
        #Open_df = pandas.read_pickle(file_path + 'Open' + file_ext)
        
        
        #%% get stooq data
        
        #data = pandas.read_pickle(r'C:\Users\srbur\stonks\stooq_data_050521.pkl')
        """
        data_volume = data.loc[:, data.columns.get_level_values(1) == 'VOLUME']
        data_volume_sum = data_volume.sum(level=[0])
        set_missing = set(Open_df.columns) - set( [i[0] for i in data_volume_sum.columns] )
        set_union = np.sort(list( set(Open_df.columns) - set_missing ))
        
        stooq = data_volume_sum[ set_union ].iloc[:-3,]
        
        
        
        #%%
        
        #Import Return Data
        Open_df = pandas.read_pickle(file_path + 'Open' + file_ext)[ set_union ].iloc[130:,] #Open from GetRetData.py #011421
        Close_df = pandas.read_pickle(file_path + 'Close' + file_ext)[ set_union ].iloc[130:,] #Close from GetRetData.py
        #Import Volume Data
        Volume_df = stooq
        #Import High Data
        High_df = pandas.read_pickle(file_path + 'High' + file_ext)[ set_union ].iloc[130:,]
        # Import Low Data
        Low_df = pandas.read_pickle(file_path + 'Low' + file_ext)[ set_union ].iloc[130:,]
        """
        
        
        
        
        
        
        """ Get the tickers and dates by index! """
        self.Tickers = Open_df.columns.values              #pandas.core.indexes.base.Index
        self.Dates = Open_df.index.values                  #pandas.core.indexes.datetimes.DatetimeIndex
        # put in a .values to get an array!
        
        """ Convert all the raw data pickles into numpy arrays! "Unassociate" the data from tickers and dates """
        self.Open = np.log(Open_df.values)                          #numpy.ndarray
        self.Close = np.log(Close_df.values)
        self.Volume = np.log(Volume_df.values)
        self.High = np.log(High_df.values)
        self.Low = np.log(Low_df.values)
        
        self.MktCap = np.log(mkt_cap_df.values)
        self.AvgVol = np.log(Avg_Volume_df.values)
        self.VolDif = np.log(Volume_df.values - volume_5min_df.values)
        # TODO MktCap is now daily now! we need to update this in the datapoints
        
        
        """ Fix Stock Splits in the Allocate_DataPoints function!
            (Remove these days from the data fed into the ML, in the datapoints array)
            GOED: day 7, stock 2510
        
        """
        
        
        
        """ Get the fundamental data for market cap! """
        '''
        Fundamental_df = pandas.read_pickle(fund_path).T
        Fundamental_df[['marketCap']] = Fundamental_df[['marketCap']].apply(pandas.to_numeric, errors='coerce')
        Fundamental_df = Fundamental_df.T
        '''
        
        """ Get the market cap into an array! """
        '''
        # Log requires float input, cast the numpy array as a float!
        self.MktCap = np.log(Fundamental_df.loc[['marketCap']].values[0].astype(float))  #numpy.ndarray
        # type(FileMama.MktCap)
        '''
    
    # Creates a new blank data day entry in all the relevant raw arrays
    # ndarray vstack and concatenate (basically same thing) take an input parameter
    # that is a list of all the arrays we want to combine!
    def New_Day(self) :
        self.Open = np.vstack((self.Open, np.array([math.nan for i in range(len(self.Tickers))])))
        self.Close = np.vstack((self.Close, np.array([math.nan for i in range(len(self.Tickers))])))
        self.Volume = np.vstack((self.Volume, np.array([math.nan for i in range(len(self.Tickers))])))
        self.High = np.vstack((self.High, np.array([math.nan for i in range(len(self.Tickers))])))
        self.Low = np.vstack((self.Low, np.array([math.nan for i in range(len(self.Tickers))])))
        
        # TODO Would need to add functionality here to handle the new data like daily market cap, avg vol, and 5 min vol!
        
        # TODO potentially add on to the dates... maybe, if we want to... or just don't
        # Even keep track of the dates at all
        
        # TODO we're actually gonna want to make everything in the FileMama
        # a pandas dataframe, or a dictionary, maybe, because it's only set up once with 
        # the raw data, then exported as a list or an array
    
    # Add temporary data so we can make a fast decision
    # The defualt values for all args will be None (maybe could also do nan, but i think none is better)
    # (The origional array is filled with nans if undefined)
    
    # Allows us to specify input in terms of a ticker, DataPointInfo object, or DataPoint array
    def Update_DataPoint(self, Ticker = None, DataPointInfo = None, DataPoint = None, Open = None, Close = None, Volume = None, High = None, Low = None) :
        # If we got a ticker input
        if not (Ticker is None) :
            # Find the index of the ticker we're adding, np where allows you to find
            # locations in an array that satisfy the function!
            tick_location = np.where(self.Tickers == Ticker)
            if len(tick_location[0]) :
                tick_index = tick_location[0][0]
                day_index = -1
            else :
                #print("Ticker '" + Ticker + "' not found")
                return None
        elif not (DataPointInfo is None) :
            tick_index = DataPointInfo.tick
            day_index = DataPointInfo.day
        elif not (DataPoint is None) :
            tick_index = round(DataPoint[10])
            day_index = round(DataPoint[0])
        else :
            print("You need to specify a datapoint array or object, or ticker")
            return None
        
        if not (Open is None) :
            self.Open[day_index][tick_index] = np.log(Open)
        if not (Close is None) :
            self.Close[day_index][tick_index] = np.log(Close)
        if not (Volume is None) :
            self.Volume[day_index][tick_index] = np.log(Volume)
        if not (High is None) :
            self.High[day_index][tick_index] = np.log(High)
        if not (Low is None) :
            self.Low[day_index][tick_index] = np.log(Low)
            
        # TODO Would need to add functionality here to handle the new data like daily market cap, avg vol, and 5 min vol!
                
        # Return the negative index to the datapoint, or stock
        #return tick_index - len(self.Tickers)
        return (day_index, tick_index)


class DataPointSource :
    # Define basic info about each of the dimensions, independent variables only!
    dim_names =  ["day", "market cap", "int-dump", "volume", "int-vol",        "prev night change", "volitility", "volume trend"]
    dim_pct_ch = [False, False,        True,       False,    True,              True,                True,         True]
    dim_exp =    [False, True,         True,       True,     True,              True,                True,         True]
    
    # Keep track of the indexes of each metric
    
    # 0 is day index, has to be first
    Day_Index = 0
    
    
    MarketCap = 1
    IntDump = 2
    Volume = 3
    VolumeAvg = 3
    PrevDay = 3
    IntVol = 4
    PrevNight = 5
    Volitility = 6
    VolumeFract = 6
    HighOpen = 6
    VolumeTrend = 7
    OvernightRetLr = 8
    OvernightRetTrans = 9
    
    
    # 10 is tick index, has to be last
    Tick_Index = 10
    
    # Ones we're not currently using, but we have functions for
    IntScore = 2
    IntGain = 2
    CloseClose = 2
    IntDrop = 4
    OvernightRet = 9
    DayCh = 6
    
    
    # Number of curve dimensions
    n_dims = len(dim_names)
    sqrt_n_dims = math.sqrt(n_dims)
    
    def __init__(self, file_mama) :
        self.HighClamp = 5
        self.file_mama = file_mama
        self.last_real_day_index = len(file_mama.Dates)-1
        
    def Allocate_DataPoints(self, start_day_index) :
        # Get the file mama
        file_mama = self.file_mama
        
        # Create a new list of the datapoint objects (used to calculate
        # the datapoint parameters we want)
        # The way this is currently set up, calling this function again
        # Will ditch the old pointer and make a new pointer, so we would
        # have to reassign references to the new pointer if we got it
        # once and then got it again
        
        # Create a new list of the datapoints
        self.DataPoints = []
        
        # Create a list allowing us to tab the index by day
        self.DayIndexes = [0 for i in range(start_day_index)]
        self.DataPointObjects = [[] for i in range(start_day_index)]
        current_index = 0
        
        # Fill out the list of datapoints
        # Cycle through the days we're dealing with
        for day_index in range(start_day_index, len(file_mama.Open)) :
            self.DayIndexes.append(current_index)
            # keep track of the current day objects
            currentDayObjects = []
            # Cycle through the tickers we're dealing with
            for tick_index in range(len(file_mama.Open[day_index])) :
                
                """ Fix Stock Splits in the Allocate_DataPoints function!
                    (Remove these days from the data fed into the ML, in the datapoints array)
                    GOED: day 7, stock 2510
                
                """
                if not ((day_index == 7) & (tick_index == 2510)) :
                    # Create a new datapoint at the day and ticker
                    DataPointObject = DataPointInfo(self, day_index, tick_index)
                    # Add it on to the list of datapoint objects
                    #self.DataPointObjects.append(DataPointObject)
                    # ADD IT ON AS A REFERENCE! This actually works by default with
                    # objects i think!
                    currentDayObjects.append(DataPointObject)
                    # Add it on to the list of datapoints
                    # ADD IT ON AS A REFERENCE!
                    # ndarray_instance.view() returns a reference
                    self.DataPoints.append(DataPointObject.DataPoint.view())
                    current_index += 1
                
            self.DataPointObjects.append(currentDayObjects)
    
    def Update_DataPoints(self) :
        for DataPointObjectRow in self.DataPointObjects :
            for DataPointObject in DataPointObjectRow :
                DataPointObject.Update_DataPoint()
    
    def Get_DataPoint_Info(self, datapoint_index_pair) :
        # We could match the datapoint with its corresponding DataPointObject
        # Or, we can just return the ticker and date for now
        #day_index = datapoint_type_list[0]
        #tick_index = datapoint_type_list[self.n_dims + self.extra_dims]
        day_index = datapoint_index_pair[0]
        tick_index = datapoint_index_pair[1]
        return self.DataPointObjects[day_index][tick_index]
    
    def Get_Ticker(self, datapoint_index_pair) :
        return self.file_mama.Tickers[datapoint_index_pair[1]]
        
        """
        ticker = self.file_mama.Tickers[tick_index]
        
        if day_index < len(self.file_mama.Dates) :
            date = self.file_mama.Dates[day_index]
        else :
            date = "New Date"
        
        return (ticker, date)
        """
    
    def Get_Point_Index(self, point) :
        day_index = round(point[self.Day_Index])
        tick_index = round(point[self.Tick_Index])
        return (day_index, tick_index)
    
    # Here's some random functions for printing that got thrown in...
    # Some should get organized and probably put somewhere else, actually
    # They do kinda belong around here since they depend on the dimension configs
    # Print out the current center and spreads
    def printCenSpr(self, cen_spr) :
        dim_names = self.dim_names
        for i in range(len(dim_names)) :
            dim_name = dim_names[i]
            print(dim_name + ": (" + Round_SF(cen_spr[i][0], 4) + ", " + Round_SF(cen_spr[i][1], 4) + ")")
            
    def printCenSprArray(self, cen_spr) :
        print("[", end = "")
        for i in range(self.n_dims) :
            if i :
                print(", ", end = "")
            print("[" + Round_SF(cen_spr[i][0], 4) + ", " + Round_SF(cen_spr[i][1], 4) + "]", end = "")
        print("]")
    
    def convertToLinRange(self, cen_spr, spr_thresh) :
        lin_range = []
        cs_index = 0
        for c_s in cen_spr :
            center = c_s[0]
            spread = c_s[1] * spr_thresh
            
            min_ran = center - spread
            max_ran = center + spread
            
            # Skip the exponential for the day index
            if self.dim_exp[cs_index] :
                min_ran = math.exp(min_ran)
                max_ran = math.exp(max_ran)
                
            if self.dim_pct_ch[cs_index] :
                min_ran -= 1
                max_ran -= 1
            
            #min_ran = Round_SF(min_ran, 4)
            #max_ran = Round_SF(max_ran, 4)
            
            lin_range.append([min_ran, max_ran])
            
            cs_index += 1
        return lin_range


class DataPointInfo :
    # Define the raw parameters used
    
    def __init__(self, data_point_source, day_index, tick_index) :
        self.dps = data_point_source
        self.fm = data_point_source.file_mama
        
        self.day = day_index
        self.tick = tick_index
        
        # Allocate float (python type for double) array
        self.DataPoint = np.ndarray(shape=(self.dps.Tick_Index+1,), dtype=float)
    
    # Define the parameters calculated
    def Int_Dump(self) :
        Close = self.fm.Close[self.day][self.tick]
        High  = self.fm.High [self.day][self.tick]
        
        self.DataPoint[self.dps.IntDump] = Close - High
    
    def Int_Vol(self) :
        High  = self.fm.High [self.day][self.tick]
        Low   = self.fm.Low  [self.day][self.tick]
        
        self.DataPoint[self.dps.IntVol] = High - Low
    
    def Int_Drop(self) :
        Low  = self.fm.Low  [self.day][self.tick]
        Open  = self.fm.Open [self.day][self.tick]
        
        self.DataPoint[self.dps.IntDrop] = Low - Open
    
    def Int_Gain(self) :
        Low  = self.fm.Low  [self.day][self.tick]
        Close  = self.fm.Close [self.day][self.tick]
        
        self.DataPoint[self.dps.IntGain] = Close - Low
    
    def High_Open(self) :
        High  = self.fm.High [self.day][self.tick]
        Open  = self.fm.Open [self.day][self.tick]
        
        self.DataPoint[self.dps.HighOpen] = High - Open
    
    def Prev_Day(self) :
        PrevOpen  = self.fm.Open [self.day-1][self.tick]
        PrevClose = self.fm.Close[self.day-1][self.tick]
        
        self.DataPoint[self.dps.PrevDay] = PrevClose - PrevOpen
    
    def Prev_Day_Vol(self) :
        PrevHigh  = self.fm.High [self.day-1][self.tick]
        PrevLow   = self.fm.Low  [self.day-1][self.tick]
        
        self.DataPoint[self.dps.PrevDay] = PrevHigh - PrevLow
    
    def Int_Score(self) :
        High  = self.fm.High [self.day][self.tick]
        Low   = self.fm.Low  [self.day][self.tick]
        Close = self.fm.Close[self.day][self.tick]
        
        self.DataPoint[self.dps.IntScore] = (High - Close) / (High - Low)
        
    def Prev_Night(self) :
        Open  = self.fm.Open [self.day][self.tick]
        PrevClose = self.fm.Close[self.day-1][self.tick]
        
        self.DataPoint[self.dps.PrevNight] = Open - PrevClose
    
    def Day_Ch(self) :
        Close = self.fm.Close[self.day][self.tick]
        Open  = self.fm.Open [self.day][self.tick]
        
        self.DataPoint[self.dps.DayCh] = (Close - Open)# * max(1, np.exp(Open))
    
    def Volitility(self) :
        # numpy.average allows weights!
        # Closes = self.fm.Close[self.day-9:, self.tick]
        # Potentially try looking at historic overnight volitility too!
        closem0 = self.fm.Close[self.day][self.tick]
        closem1 = self.fm.Close[self.day-1][self.tick]
        closem2 = self.fm.Close[self.day-2][self.tick]
        closem3 = self.fm.Close[self.day-3][self.tick]
        closem4 = self.fm.Close[self.day-4][self.tick]
        closem5 = self.fm.Close[self.day-5][self.tick]
        closem6 = self.fm.Close[self.day-6][self.tick]
        closem7 = self.fm.Close[self.day-7][self.tick]
        closem8 = self.fm.Close[self.day-8][self.tick]
        closem9 = self.fm.Close[self.day-9][self.tick]
        closem10 = self.fm.Close[self.day-10][self.tick]
        
        self.DataPoint[self.dps.Volitility] = np.sqrt( abs(closem0 - closem1)     
                                                + abs(closem1 - closem2)
                                                + abs(closem2 - closem3)
                                                + abs(closem3 - closem4)     
                                                + abs(closem4 - closem5)
                                                + abs(closem5 - closem6)
                                                + abs(closem6 - closem7)
                                                + abs(closem7 - closem8)
                                                + abs(closem8 - closem9)
                                                + abs(closem9 - closem10) ) / 10
    
    def Volitility2(self) :
        # numpy.average allows weights!
        # Closes = self.fm.Close[self.day-9:, self.tick]
        # Potentially try looking at historic overnight volitility too!
        closem0 = self.fm.Close[self.day][self.tick]
        closem1 = self.fm.Close[self.day-1][self.tick]
        closem2 = self.fm.Close[self.day-2][self.tick]
        closem3 = self.fm.Close[self.day-3][self.tick]
        closem4 = self.fm.Close[self.day-4][self.tick]
        closem5 = self.fm.Close[self.day-5][self.tick]
        
        self.DataPoint[self.dps.Volitility] = np.sqrt( abs(closem0 - closem1)     
                                                + abs(closem1 - closem2)
                                                + abs(closem2 - closem3)
                                                + abs(closem3 - closem4)     
                                                + abs(closem4 - closem5) ) / 5
    
    def Volitility3(self) :
        # numpy.average allows weights!
        # Closes = self.fm.Close[self.day-9:, self.tick]
        # Potentially try looking at historic overnight volitility too!
        #closem0 = self.fm.Close[self.day][self.tick]
        closem1 = self.fm.Close[self.day-1][self.tick]
        closem2 = self.fm.Close[self.day-2][self.tick]
        closem3 = self.fm.Close[self.day-3][self.tick]
        closem4 = self.fm.Close[self.day-4][self.tick]
        closem5 = self.fm.Close[self.day-5][self.tick]
        closem6 = self.fm.Close[self.day-6][self.tick]
        closem7 = self.fm.Close[self.day-7][self.tick]
        closem8 = self.fm.Close[self.day-8][self.tick]
        closem9 = self.fm.Close[self.day-9][self.tick]
        closem10 = self.fm.Close[self.day-10][self.tick]
        
        self.DataPoint[self.dps.Volitility] = np.sqrt( 0     
                                                + abs(closem1 - closem2)
                                                + abs(closem2 - closem3)
                                                + abs(closem3 - closem4)     
                                                + abs(closem4 - closem5)
                                                + abs(closem5 - closem6)
                                                + abs(closem6 - closem7)
                                                + abs(closem7 - closem8)
                                                + abs(closem8 - closem9)
                                                + abs(closem9 - closem10) ) / 9
    
    def VolumeFract(self) :
        vol_total = self.fm.Volume[self.day][self.tick]
        vol_diff  = self.fm.VolDif[self.day][self.tick]
        
        self.DataPoint[self.dps.VolumeFract] = vol_diff - vol_total
    
    def VolumeFractTrend(self) :
        volm0 = self.fm.VolDif[self.day][self.tick]
        volm1 = self.fm.VolDif[self.day-1][self.tick]
        
        self.DataPoint[self.dps.VolumeFract] = volm0 - volm1
    
    def Close_Close(self) :
        Closem0 = self.fm.Close[self.day][self.tick]
        Closem1  = self.fm.Close [self.day-1][self.tick]
        
        self.DataPoint[self.dps.CloseClose] = Closem0 - Closem1
    
    def Volume(self) :
        Volume = self.fm.Volume[self.day][self.tick]
        
        self.DataPoint[self.dps.Volume] = Volume
    
    def Volume_Avg(self) :
        volm0 = self.fm.Volume[self.day][self.tick]
        volm1 = self.fm.Volume[self.day-1][self.tick]
        volm2 = self.fm.Volume[self.day-2][self.tick]
        volm3 = self.fm.Volume[self.day-3][self.tick]
        volm4 = self.fm.Volume[self.day-4][self.tick]
        volm5 = self.fm.Volume[self.day-5][self.tick]
        
        self.DataPoint[self.dps.VolumeAvg] = (volm0 + volm1 + volm2 + volm3 + volm4 + volm5) / 6
        
    
    def Volume_Trend(self) :
        # TODO potentially figure out a better way of doing this
        volm0 = self.fm.Volume[self.day][self.tick]
        volm1 = self.fm.Volume[self.day-1][self.tick]
        volm2 = self.fm.Volume[self.day-2][self.tick]
        volm3 = self.fm.Volume[self.day-3][self.tick]
        
        self.DataPoint[self.dps.VolumeTrend] = ( 1.00 * (volm0 - volm1)     
                                               + 0.61 * (volm1 - volm2)
                                               + 0.37 * (volm2 - volm3) )
    
    # Must be calculated after volume average
    def Volume_Trend2(self) :
        volm0 = volm0 = self.fm.Volume[self.day][self.tick]
        volume_avg = self.DataPoint[self.dps.VolumeAvg]
        
        self.DataPoint[self.dps.VolumeTrend] = volm0 - volume_avg
    
    def Volume_Trend3(self) :
        volm0 = self.fm.Volume[self.day][self.tick]
        volm1 = self.fm.Volume[self.day-1][self.tick]
        volm2 = self.fm.Volume[self.day-2][self.tick]
        volm3 = self.fm.Volume[self.day-3][self.tick]
        volm4 = self.fm.Volume[self.day-4][self.tick]
        volm5 = self.fm.Volume[self.day-5][self.tick]
        
        self.DataPoint[self.dps.VolumeTrend] = ( 1.00 * (volm0 - volm1)     
                                               + 0.61 * (volm1 - volm2)
                                               + 0.37 * (volm2 - volm3)
                                               + 0.22 * (volm3 - volm4)
                                               + 0.14 * (volm4 - volm5))
    
    def Volume_Trend3b(self) :
        volm0 = self.fm.Volume[self.day][self.tick]
        volm1 = self.fm.Volume[self.day-1][self.tick]
        volm2 = self.fm.Volume[self.day-2][self.tick]
        volm3 = self.fm.Volume[self.day-3][self.tick]
        volm4 = self.fm.Volume[self.day-4][self.tick]
        volm5 = self.fm.Volume[self.day-5][self.tick]
        
        vol_trend = ( 1.00 * (volm0 - volm1)     
                    + 0.61 * (volm1 - volm2)
                    + 0.37 * (volm2 - volm3)
                    + 0.22 * (volm3 - volm4)
                    + 0.14 * (volm4 - volm5))
        
        if vol_trend > 0 :
            self.DataPoint[self.dps.VolumeTrend] = np.log(1+vol_trend)
        else :
            self.DataPoint[self.dps.VolumeTrend] = vol_trend
    
    def Volume_Trend4(self) :
        volm0 = self.fm.Volume[self.day][self.tick]
        volm1 = self.fm.Volume[self.day-1][self.tick]
        volm2 = self.fm.Volume[self.day-2][self.tick]
        volm3 = self.fm.Volume[self.day-3][self.tick]
        volm4 = self.fm.Volume[self.day-4][self.tick]
        volm5 = self.fm.Volume[self.day-5][self.tick]
        volm6 = self.fm.Volume[self.day-6][self.tick]
        volm7 = self.fm.Volume[self.day-7][self.tick]
        
        self.DataPoint[self.dps.VolumeTrend] = ( 1.00 * (volm0 - volm1)     
                                               + 0.61 * (volm1 - volm2)
                                               + 0.37 * (volm2 - volm3)
                                               + 0.22 * (volm3 - volm4)
                                               + 0.14 * (volm4 - volm5)
                                               + 0.08 * (volm5 - volm6)
                                               + 0.05 * (volm6 - volm7))
    
    def Volume_Trend_Avg(self) :
        volm0   = self.fm.Volume[self.day][self.tick]
        volavg  = self.fm.AvgVol[self.day][self.tick]
        
        self.DataPoint[self.dps.VolumeTrend] = volm0 - volavg
    
    def Market_Cap(self) :
        MarketCap = self.fm.MktCap[self.tick]
        
        self.DataPoint[self.dps.MarketCap] = MarketCap
    
    def Market_Cap_Daily(self) :
        MarketCap = self.fm.MktCap[self.day][self.tick]
        
        self.DataPoint[self.dps.MarketCap] = MarketCap
    
    def Overnight_Ret(self) :
        # Only calculate this one if we're not on the last day! Otherwise nan
        if self.day < self.dps.last_real_day_index :
            NextOpen  = self.fm.Open [self.day+1][self.tick]
            Close = self.fm.Close[self.day][self.tick]
            
            self.DataPoint[self.dps.OvernightRet] = NextOpen - Close
        else :
            self.DataPoint[self.dps.OvernightRet] = math.nan
    
    def Overnight_Ret_Trans(self) :
        # Only calculate this one if we're not on the last day! Otherwise nan
        if self.day < self.dps.last_real_day_index :
            NextOpen  = self.fm.Open [self.day+1][self.tick]
            Close = self.fm.Close[self.day][self.tick]
            
            # Factor in 1% total transaction loss
            self.DataPoint[self.dps.OvernightRetTrans] = NextOpen - Close - 0.008#0.01
        else :
            self.DataPoint[self.dps.OvernightRetTrans] = math.nan
    
    # Must calculate after Overnight_Ret!
    def Overnight_Ret_Lr(self) :
        OvernightRet = self.DataPoint[self.dps.OvernightRet]
        
        self.DataPoint[self.dps.OvernightRetLr] = OvernightRet
        
        if OvernightRet > 0 :
            high_clamp = self.dps.HighClamp
            
            self.DataPoint[self.dps.OvernightRetLr] = np.log(1+OvernightRet*high_clamp)/high_clamp
    
    
    def Update_DataPoint(self) :
        # Calculate all the calculated parameters we want to calculate!
        #self.Int_Dump()
        #self.Close_Close()
        self.Int_Gain()
        #self.Int_Score()
        self.Int_Vol()
        #self.Int_Drop()
        
        self.Prev_Night()
        #self.Day_Ch()
        #self.Volitility()
        #self.Volitility2()
        #self.VolumeFract()
        #self.Day_Ch()
        self.High_Open()
        #self.VolumeFractTrend()
        #self.Volitility3()
        #self.Volume()
        self.Prev_Day()
        #self.Prev_Day_Vol()
        #self.Volume_Avg()
        #self.Volume_Trend()
        #self.Volume_Trend2()
        #self.Volume_Trend3()
        self.Volume_Trend_Avg()
        #self.Volume_Trend3b()
        #self.Volume_Trend4()
        #self.Market_Cap()
        self.Market_Cap_Daily()
        
        self.Overnight_Ret()
        #self.Overnight_Ret_Trans()
        self.Overnight_Ret_Lr()
        
        # This covers all the derived parameters.
        # We also need to set the fixed parameters like day and stock
        self.DataPoint[self.dps.Day_Index] = float(self.day)
        self.DataPoint[self.dps.Tick_Index] = float(self.tick)
        
        # Put all the desired parameters into the datapoint in the order we want!
        # This means that all of these things are stored twice in memory,
        # But it seems worth it for the code organization
        
        # TODO: Potentially have the last index point back to the 
        # DataPointInfo Object

"""
FileMama = FileMama(-5)

FileMama.New_Day()

FileMama.Update_DataPoint('A', math.nan, 2, 1000, 2.2, 0.9)

DataPointSource = DataPointSource(FileMama)

DataPointSource.Allocate_DataPoints(3)

print("DataPointSource.DataPointObjects[3][0].day: " + str(DataPointSource.DataPointObjects[3][0].day))

DataPointSource.Update_DataPoints()
"""
#DataPointSource.DataPoints[0]

#DataPointSource.Get_DataPoint_Info(DataPointSource.DataPoints[0])

#DataPointSource.DataPoints[DataPointSource.DayIndexes[5]]

#DataPointSource.Get_DataPoint_Info(DataPointSource.DataPoints[DataPointSource.DayIndexes[5]])

#DataPointSource.DataPointObjects[3][0] -> Get datapointobject from day and ticker

#DataPointSource.DataPointObjects[3][0].DataPoint


    
