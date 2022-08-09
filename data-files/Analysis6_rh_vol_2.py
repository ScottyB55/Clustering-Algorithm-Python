# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 10:12:10 2021

@author: twbur
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 15:50:50 2021

@author: twbur
"""

#file_path = r'C:\Users\twbur\Documents\CodePython\Robinhood\rh_volume\rh_volume_data.pkl'
#file_path_5min = r'C:\Users\twbur\Documents\CodePython\Robinhood\rh_5min\rh_5min_data.pkl'
file_path = r'rh_volume_data.pkl'
file_path_5min = r'rh_5min_data.pkl'

#Purpose: analyze relationship between RH positioning and equity returns. Use pct. return, price, volume as
#proxies for popular RH stocks.

# Cleaned, better commented version of Analysis1 code

# Import packages and data

import numpy as np
import pandas
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
data = pandas.read_pickle(file_path)
data_5min = pandas.read_pickle(file_path_5min)

#%% get volume data
data_volume = data_5min.loc[:, data_5min.columns.get_level_values(1) == 'volume']
data_volume_sum = data_volume.sum(level=[0])
set_missing = set( [i[0] for i in data.columns] ) - set( [i[0] for i in data_volume_sum.columns] )
set_union = np.sort(list( set([i[0] for i in data.columns]) - set_missing ))

volume_5min = data_volume_sum[ set_union ]

#%%

#Import Return Data
Open_df = data.loc[:, data.columns.get_level_values(1) == 'open'][ set_union ]
Close_df_0 = data.loc[:, data.columns.get_level_values(1) == 'market_cap'][ set_union ].values / data.loc[:, data.columns.get_level_values(1) == 'shares_outstanding'][ set_union ].values
Close_df = pandas.DataFrame(Close_df_0)
Open = Open_df.values
Close = Close_df.values
#Import Volume Data
Volume_df = data.loc[:, data.columns.get_level_values(1) == 'volume'][ set_union ]
Volume = Volume_df.values
#Import Avg Volume Data
Avg_Volume_df = data.loc[:, data.columns.get_level_values(1) == 'average_volume'][ set_union ]
Avg_Volume = Avg_Volume_df.values
#Import High Data
High_df = data.loc[:, data.columns.get_level_values(1) == 'high'][ set_union ]
High = High_df.values
# Import Low Data
Low_df = data.loc[:, data.columns.get_level_values(1) == 'low'][ set_union ]
Low = Low_df.values

# Import Fundamental Data
mkt_cap_df = data.loc[:, data.columns.get_level_values(1) == 'market_cap'][ set_union ]
mkt_cap = mkt_cap_df.values

#%% Get data that tries to predict high RH volume

#rolling avg volume pct chg
Volume_pct_ch = Volume/Avg_Volume-1

#volume pct chg
# Volume_pct_ch_df = Volume_df.pct_change()
# Volume_pct_ch = Volume_pct_ch_df.values

#volume relative to 5min volume pct chg
#Volume_alt = Volume_df.values/volume_5min.values-1
Volume_alt = (Volume_df.values-volume_5min.values)/Avg_Volume


Close_pct_ch_df = Close_df.pct_change()
Close_pct_ch = Close_pct_ch_df.values

pr_diff = Close-Open
pr_pct_ch = pr_diff/Open

high_close_diff = Close - High
high_close_pct_ch = high_close_diff/High

high_low_diff = Low - High
high_low_pct_ch = high_low_diff/High

night_pct_chg = Close_pct_ch - pr_pct_ch


#%% NEW TEST NEW TEST NEW TEST NEW TEST
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#Test only using market data (no rh data b/c robinhood data is gone now)

#%% Analyze the effects of volume and daily price change on overnight price chg

# new ~~~~~~~~~~~~~~~~~~~~~~~~~~
#Conditions to be met - lol
pos_bool3 = ( Close_pct_ch > -.07 ) & ( Close_pct_ch < .07 ) #close pct change high enough; low=bad
pos_bool4 = ( pr_pct_ch > -.07 ) & ( pr_pct_ch < .07 ) #must have actually gone up intraday
pos_bool5 = ( Volume_alt > 8 ) & ( Volume_alt < 400 ) #change in volume cannot be too low or high (try to screen out real events)
pos_bool6 = ( mkt_cap < .75*1000000000 ) & ( Close < 20 ) #penny stock screen. mkt_cap in billions
pos_bool7 = Volume > .8*1000000 #liquidity requirement. in millions
pos_bool8 = high_close_pct_ch > -.09 #stock did not massively dump intraday
pos_bool12 = (high_low_pct_ch > -.1)
pos_bool13 = ( night_pct_chg > -.09 ) & ( night_pct_chg < .09 )
#Combine conditions
pos_bool = pos_bool3*pos_bool4*pos_bool5*pos_bool6*pos_bool7*pos_bool8
pos_bool = pos_bool[:-1]

# old ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Conditions to be met
# pos_bool3 = ( Close_pct_ch > .14 ) & ( Close_pct_ch < 1 ) #close pct change high enough; low=bad
# pos_bool4 = pr_pct_ch > .09 #must have actually gone up intraday
# pos_bool5 = ( Volume_pct_ch > .5 ) & ( Volume_pct_ch < 10 ) #change in volume cannot be too low or high (try to screen out real events)
# pos_bool6 = (Close < 12) #( mkt_cap < 2*1000000000 ) #& (Close > 2) #penny stock screen
# pos_bool7 = Volume > 110 #liquidity requirement
# pos_bool8 = high_close_pct_ch > -.5 #stock did not massively dump intraday
# #Combine conditions
# pos_bool = pos_bool3*pos_bool4*pos_bool5*pos_bool6*pos_bool7*pos_bool8
# pos_bool = pos_bool[:-1]

#Strategy
hp_pr_diff = Open[1:,:] - Close[:-1,:] #or: Close[1:,:] - Close[:-1,:] #
hp_pr_pct_ch = hp_pr_diff/Close[:-1,:]
strat_rets = hp_pr_pct_ch*pos_bool #Close_pct_ch[1:,:]*pos_bool gives close-to-close
strat_rets = np.where(np.isnan(strat_rets), 0, strat_rets)

#print('RESULT: ' + str( sum(sum(strat_rets))/sum(sum(pos_bool)) ) + ', in: ' + str( sum(sum(pos_bool))/len(strat_rets)) + ' trades/day' )


#Strategy returns
L = 1 #fraction of account dedicated to trading stock (leverage)
strat_optimal_compound_rets = np.sum(strat_rets,axis=1)/np.sum(pos_bool,axis=1) - 0.0 #Creates daily portfolio return, assuming equal weights
daily_rets = strat_optimal_compound_rets #keep for future analysis
strat_optimal_compound_rets = np.where(np.isnan(strat_optimal_compound_rets), 0, strat_optimal_compound_rets)*L+1
strat_optimal_compound_rets = np.cumprod(strat_optimal_compound_rets)

#plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_yscale('linear')
#ax.plot(Open_df.index[1:][:], strat_optimal_compound_rets[:])
day_start = -25
ax.plot(Open_df.index[1:][:], strat_optimal_compound_rets[:])
plt.xticks(rotation=45)
plt.title('growth of $1')

print('RESULT: ' + str( np.sum(np.sum(strat_rets))/np.sum(np.sum(pos_bool)) ) + ', in: ' + str( np.sum(np.sum(pos_bool))/len(strat_rets)) + ' trades/day' )
print('Compound Ret Result: ' + str( strat_optimal_compound_rets[-1] ))



#%% AUDIT TRADES - see specific trades the strategy recommends and the results of those trades

#See specific trades
xy = np.nonzero(pos_bool == 1)
date_arr = []
tic_arr = []
ret_arr = []
day_ret_arr = []
volume_ch_arr = []
volume_alt_arr = []
intraday_ret_arr = []

for x in xy[0]:
    date_arr.append(Open_df.index[x])
for y in xy[1]:
    tic_arr.append(Open_df.columns[y][0])
for i in np.arange(len(date_arr)):
    ret_arr.append(strat_rets[ xy[0][i] , xy[1][i] ])
    volume_ch_arr.append(Volume_pct_ch[ xy[0][i] , xy[1][i] ])
    volume_alt_arr.append(Volume_alt[ xy[0][i] , xy[1][i] ])
    day_ret_arr.append(Close_pct_ch[ xy[0][i] , xy[1][i] ])
    intraday_ret_arr.append(pr_pct_ch[ xy[0][i] , xy[1][i] ])

trades = pandas.DataFrame(np.transpose([date_arr,tic_arr,ret_arr,day_ret_arr,intraday_ret_arr,volume_ch_arr,volume_alt_arr]))
trades.columns = ["Buy Date","Tic","Overnight Ret","day_ret","intraday_ret","volume_ch","volume_alt"]
#trades.to_csv('strategy_recommended_trades.csv')
trades