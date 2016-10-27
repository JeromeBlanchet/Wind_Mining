import datetime
import weatherMethods


# UTC Flight Departure Time
departureTime = datetime.datetime(2013, 2, 3, 12, 0)

print(departureTime)

# Flight Distance (km)
flightdistance = 500

# This function downloads all the data needed for the flight
meteoDT, winds, lvls, timeTags = weatherMethods.getWeather(departureTime, flightdistance)

#Current position
lon = -90
lat = 43

# Current altitude (ft)
alt = 30000

# UTC time at currentpostion
t = datetime.datetime(2013, 2, 3, 16, 30)

# Last position
lonl = -89.9
latl= 43.1

# This function gives wind value at current position (km/h)
wind  = weatherMethods.getWind(lon, lat, alt, t, meteoDT, lonl, latl, timeTags, winds, lvls)

print(wind)

