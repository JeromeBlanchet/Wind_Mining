
# coding: utf-8

# In[133]:

import datetime
import weatherMethods
import timeMethods
import os
import numpy as np
import pandas as pd
import time
from tools import press, proxilvl, g
import pickle
from scipy.spatial import cKDTree
import sys
import math

"""
Required Files:
Valid FLight Tracks
Labeled Data
MNL_DATA
"""

DEP = 'JFK'
ARR = 'LAX'
Year = 2013
print(DEP,'-->',ARR)


# # Initialize Grid System & Construct lvl-specific KDTree
# 
# All US mainland

# In[2]:

departureTime = datetime.datetime(2013, 2, 1, 1, 30)
arrivalTime = datetime.datetime(2013, 2, 1, 6, 30)
WindClass = weatherMethods.GetWindSpeed(departureTime, arrivalTime)
Winds = WindClass.winds


# In[4]:

# Tree_Dict = {}
# i = 0
# for levels in WindClass.lvls.keys():
#     i += 1
#     Grid_lat, Grid_lon = Winds[0][0][WindClass.lvls[levels]].data(22, 52, -130, -60)[1:3]
# #     Grid_lat, Grid_lon = Winds[0][0][WindClass.lvls[levels]].data(28, 45, -99, -68)[1:3]
# #     print(Grid.shape) 145999 * 2
#     Grid = np.dstack((Grid_lat, Grid_lon))[0]
#     Tree_Dict[levels] = cKDTree(Grid)
# pickle.dump(Tree_Dict, open('Grid_KDTree.p','wb'))


# In[5]:

Tree_Dict = pickle.load(open('Grid_KDTree.p','rb'))


# # Map Trajectories

VTrackPath = os.getcwd() + '/TFMS_NEW/New_' + DEP + ARR + str(Year) + '.csv'
VTrack = pd.read_csv(VTrackPath, parse_dates=[6])
LabelData = pd.read_csv(os.getcwd() + '/TFMS_NEW/Label_' + DEP+'_' + ARR+ '_' + str(Year) + '.csv', parse_dates=[6])
CenterTraj = VTrack[VTrack.FID.isin(LabelData[LabelData.MedianID != -2].FID.values)].reset_index(drop = 1)


# In[7]:

CenterTraj['levels'] = CenterTraj['Alt'].apply(lambda x: proxilvl(x*100, WindClass.lvls))
CenterTraj['QueryIdx'] = 0
CenterTraj['QueryIdx'] = CenterTraj['QueryIdx'].astype(int)
for lvl, gp in CenterTraj.groupby('levels'):
    CenterTraj.loc[gp.index, 'QueryIdx'] = Tree_Dict[lvl].query(gp[['Lat','Lon']])[1]


# In[27]:

Tree_index = CenterTraj[['levels','FID','Lat','Lon','GroundSpeed','DT','Dist','QueryIdx']].set_index(['levels'])


# In[28]:

print('----------------Start Mapping Wind with Trajectories----------------')
MissingFID = []
st = time.time()
for i in range(LabelData.shape[0]):
    if i % 500 == 0:
        print(i, time.time() - st)
    try:
        if i == 0:
            temp_wind = Tree_index.copy()
        else:
            pass
        departureTime = LabelData.loc[i, 'Elap_Time']
        arrivalTime = departureTime + datetime.timedelta(hours = 12)
        WindClass = weatherMethods.GetWindSpeed(departureTime, arrivalTime)
        Winds = WindClass.winds
        
        u_col_name = 'u_wind_' + str(LabelData.loc[i, 'FID'])
        v_col_name = 'v_wind_' + str(LabelData.loc[i, 'FID'])
        
        temp_wind[u_col_name] = 0.0
        temp_wind[v_col_name] = 0.0
        for lvl in Tree_index.index.unique():
            try:
                u_wind = Winds[0][0][WindClass.lvls[lvl]].data(22, 52, -130, -60)[0][Tree_index.loc[lvl,'QueryIdx']]
                v_wind = Winds[0][1][WindClass.lvls[lvl]].data(22, 52, -130, -60)[0][Tree_index.loc[lvl,'QueryIdx']]
            except KeyError:
                u_wind = np.nan
                v_wind = np.nan
                
            temp_wind.set_value(lvl, u_col_name, u_wind)
            temp_wind.set_value(lvl, v_col_name, v_wind)
    except KeyboardInterrupt:
        break
    except:
        MissingFID.append(LabelData.loc[i,'FID'])
        print(sys.exc_info()[0], LabelData.loc[i,'FID'])
print(time.time() - st)
pickle.dump(MissingFID, open('MissingFID' + DEP+ARR+str(Year) + '.p','wb'))

# In[40]:

Final_wind = temp_wind.copy()
Final_wind = Final_wind.reset_index(drop = False)
Final_wind.to_csv(os.getcwd() + '/Wind_Mapping_result/' + DEP + '_'+ ARR + '_' + str(2013) + '.csv', index = False)


# In[325]:

def GetAzimuth(Final_wind1):
    Final_wind = Final_wind1.copy()
    Final_wind['azimuth'] = 0.0
    FID = []
    for rowid, row in Final_wind.iterrows():
        if int(row.FID) not in FID:
            FID.append(int(row.FID))
            latl = row.Lat
            lonl = row.Lon
        else:
            Final_wind.loc[rowid,'azimuth'] = g.inv(lonl, latl, row.Lon, row.Lat)[0]
            latl = row.Lat
            lonl = row.Lon
    return Final_wind


# In[329]:

def GetHeadWind(Final_wind):
    Final_Head_Wind = Final_wind[['FID', 'Lat', 'Lon', 'levels', 'QueryIdx',
                              'azimuth','GroundSpeed', 'DT', 'Dist']].copy()
    # Indices need to be justified
    Columns = Final_wind.columns[9:]
    for i in range(Columns.shape[0]//2):
        ColName = 'Headwind_' + Columns[i*2][7:]
        ColName1 = 'WindDist_' + Columns[i*2][7:]
        Final_Head_Wind[ColName] = (Final_wind[Columns[2*i]] * Final_wind['azimuth'].apply(lambda x: math.sin(x*math.pi/180)) + 
                                    Final_wind[Columns[2*i+1]] * Final_wind['azimuth'].apply(lambda x: math.cos(x*math.pi/180)))
        # m/s
        Final_Head_Wind[ColName1] = Final_Head_Wind[ColName] * Final_Head_Wind['DT']
    return Final_Head_Wind


# In[202]:

Final_wind = GetAzimuth(Final_wind)
Final_Head_Wind = GetHeadWind(Final_wind)

Final_Head_Wind.to_csv(os.getcwd() + '/Wind_Mapping_result/' + DEP + '_'+ ARR + '_' + str(2013) + '_WindSpeed.csv', index = False)
Final_wind.to_csv(os.getcwd() + '/Wind_Mapping_result/' + DEP + '_'+ ARR + '_' + str(2013) + '.csv', index = False)


# ### Compute Wind-related Metric

# In[319]:

IAH_BOS_MNL = pd.read_csv(os.getcwd() + '/MNL_DATA/MNL_' + DEP+'_' + ARR+ '_' + str(Year) + '_Mean.csv')


# In[287]:

MeanWindSP = Final_Head_Wind.groupby('FID')[Final_Head_Wind.columns[9::2]].mean().unstack().reset_index() # m/s
WindDist = (Final_Head_Wind.groupby('FID')[Final_Head_Wind.columns[10::2]].sum()*0.0005399568034555).unstack().reset_index() # nmi
MeanGSP = Final_Head_Wind.groupby('FID').agg({'GroundSpeed':np.mean, 'Dist': np.sum}).unstack().reset_index() # knot
# Wide to Long
WindDist['FID_Member'] = WindDist['level_0'].apply(lambda x: int(x[9:]))
WindDist.columns = ['old','FID_x', 'Wind_Dist','FID_Member']
MeanWindSP['FID_Member'] = MeanWindSP['level_0'].apply(lambda x: int(x[9:]))
MeanWindSP.columns = ['old','FID_x', 'MeanWindSpeed','FID_Member']


# In[306]:

IAH_BOS_MNL_NEW = IAH_BOS_MNL.merge(WindDist[['FID_x','FID_Member','Wind_Dist']], left_on=['FID_x','FID_Member'], right_on = ['FID_x','FID_Member'], how='left')
IAH_BOS_MNL_NEW = IAH_BOS_MNL_NEW.merge(MeanWindSP[['FID_x','FID_Member','MeanWindSpeed']], left_on=['FID_x','FID_Member'], right_on = ['FID_x','FID_Member'], how='left')


# In[316]:

IAH_BOS_MNL_NEW.to_csv(os.getcwd() + '/MNL_DATA/NEW_MNL_' + DEP+ARR+str(Year) +'.csv', index = False)

