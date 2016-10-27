# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 17:39:24 2015

@author: Emmanuel Boidot
"""
import datetime
import timeMethods
import calendar

import traceback
import sys
import os
import subprocess as commands

import pygrib

VERBOSE = 0

DATADIR = 'data'
WINDDATADIR='winddata'
RADARDATADIR='raddata'

path = '/media/liuyulin101/Elements/Wind_Data'

WEATHERDATADIR = path + '/DATA/weather_data'
FILTEREDWEATHERDATADIR = path + '/DATA/filtered_weather_data'
RAPDATADIR=WEATHERDATADIR+'/rap'
NAMDATADIR=WEATHERDATADIR+'/nam'
NAMANLDATADIR=WEATHERDATADIR+'/namanl'
GFSDATADIR=WEATHERDATADIR+'/gfs'

def get_fname(dsrc,depTime,mforecast):
    """
    Returns the name of the grib file on the NCDC server at the time requested
    
    Args:
        * 'dsrc' (string): type of source files for weather: either namanl or rap
        * 'DT (string): date and time of the grib file requested
        * 'mforecast' (string): time horizon of the weather forecast requested
    
    Outputs:
        * 'fname' (string): name of the requested grib file on the NCDC server
    """
    mdate = depTime[0:8]
    mtime = depTime[8:12]
    fname = dsrc + '/'
    fname = fname + ('rap_130_' if dsrc=='rap' else ('nam_218_' if dsrc=='nam' else 'namanl_218_'))
    fname = fname + mdate +'_'+mtime+'_'+mforecast
    fname = fname + ('.grb2' if dsrc=='rap' else '.grb')
    return fname

def download_and_filter_data(depTime,forecast='000',filetype='rap'):
    """
    Downloads the grib file from the NCDC server on the local disk
    
    Args:
        * 'DT' (string): date and time of the last weather bulletin available
        * 'forecast' (string): time horizon of the weather forecast requested
        * 'filetype' (string): type of source files for weather: either namanl or rap
        * 'verbose' (boolean): if True, print some directory paths
    """
    if not filetype in ['rap','nam','namanl']:
        raise NameError("type should be either 'rap' or 'nam'")
        return
    
    date = depTime[0:8]
    time = depTime[8:12]    
    
    ufile = DATADIR+'/wdata_'+date+'_'+time+'_u.grb'+('2' if filetype=='rap' else '')
    vfile = DATADIR+'/wdata_'+date+'_'+time+'_v.grb'+('2' if filetype=='rap' else '')
    rfile = DATADIR+'/rdata_'+date+'_'+time+'_tisrgrd.grb'+('2' if filetype=='rap' else '')

    #if not (os.path.exists(ufile) and os.path.exists(vfile) and os.path.exists(rfile)):
    if filetype in ['rap', 'nam']:
        fname = get_fname(filetype,depTime,forecast)[4:]
    else: 
        fname = get_fname(filetype,depTime,forecast)[7:]

    if not os.path.exists(FILTEREDWEATHERDATADIR+'/'+filetype+'/'+fname):
        if not os.path.exists(WEATHERDATADIR+'/'+filetype+'/'+fname):
            print('Downloading file %s ...'%fname)
            if filetype=='rap':
                mycmd = ('wget -P '+RAPDATADIR+' ftp://nomads.ncdc.noaa.gov/RUC/13km/%s/%s/'+fname)%(date[:-2],date)
            elif filetype=='namanl':
                mycmd = ('wget -P '+NAMANLDATADIR+' ftp://nomads.ncdc.noaa.gov/NAM/analysis_only/%s/%s/'+fname)%(date[:-2],date)
            elif filetype=='nam':
                mycmd = ('wget -P '+NAMDATADIR+' http://nomads.ncdc.noaa.gov/data/meso-eta-hi/%s/%s/'+fname)%(date[:-2],date)
            # only for last few days
            elif filetype=='gfs': 
                mycmd = ('wget -P '+GFSDATADIR+' ftp://nomads.ncep.noaa.gov:9090/dods/gfs/gfs%s/gfs%s_%sz')%(mdate,mdate,mforecast[1:])
            commands.getstatusoutput(mycmd)
        else:
            print('I have the file %s'%fname)
        mycmd = 'grib_filter '+'weather/filter_for_wind_and_radar_grouped_2 '+WEATHERDATADIR+'/'+filetype+'/'+fname+';'
        mycmd = mycmd+'mv '+FILTEREDWEATHERDATADIR+'/temp.grb'+('2' if filetype=='rap' else '')+' '+FILTEREDWEATHERDATADIR+'/'+filetype+'/'+fname
        commands.getstatusoutput(mycmd)
    else:
        print('I already have the filtered file %s'%fname)
        pass

    return os.path.join(os.getcwd(), ufile), os.path.join(os.getcwd(), vfile), os.path.join(os.getcwd(), rfile)



timeTags = ['000','001','002','003','006']
for month in [9]:
	for day in range(1,calendar.monthrange(2013,month)[1]+1):
		for hour in [0,6,12,18]:
			departureTime = datetime.datetime(2013,month,day,hour,0)
			depTime = timeMethods.getLastWeatherDate(departureTime, 'namanl')
			for mforecast in timeTags:
				download_and_filter_data(timeMethods.dateFromDT(depTime), mforecast, 'namanl')




