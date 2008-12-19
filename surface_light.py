import math, matplotlib

"""
Functions for gregorian calendar and estimate of surface light converted to
python on 10.06.2008 by Trond Kristiansen
Trond.Kristiansen@imr.no
"""

def s2hms(secs):
    """
    S2HMS converts seconds to integer hour,minute,seconds
    """
    sec=round(secs)
    hour=math.floor(sec/3600)
    min=math.floor((sec % 3600)/60)
    sec=round((sec % 60))

    return hour, min, sec
    
def gregorian(julian):
    
    """
    GREGORIAN  Converts digital Julian days to Gregorian calendar dates.  
           Although formally, 
           Julian days start and end at noon, here Julian days
           start and end at midnight for simplicity.
         
           In this convention, Julian day 2440000 begins at 
           0000 hours, May 23, 1968.
    
         Usage: [gtime]=gregorian(julian) 
    
            julian... input decimal Julian day number
    
            gtime is a six component Gregorian time vector
              i.e.   gtime=[yyyy mo da hr mi sec]
                     gtime=[1989 12  6  7 23 23.356]
    
            yr........ year (e.g., 1979)
            mo........ month (1-12)
            d........ corresponding Gregorian day (1-31)
            h........ decimal hours
    
      calls S2HMS
     
        julian=julian+5.e-9;    % kludge to prevent roundoff error on seconds
    
          if you want Julian Days to start at noon...
          h=rem(julian,1)*24+12;
          i=(h >= 24);
          julian(i)=julian(i)+1;
          h(i)=h(i)-24;
    """
    secs=(julian % 1)*24*3600

    j = math.floor(julian) - 1721119
    inn = 4*j -1
    y = math.floor(inn/146097)
    j = inn - 146097*y
    inn = math.floor(j/4)
    inn = 4*inn +3
    j = math.floor(inn/1461);
    d = math.floor(((inn - 1461*j) +4)/4)
    inn = 5*d -3
    m = math.floor(inn/153)
    d = math.floor(((inn - 153*m) +5)/5)
    y = y*100 +j
    mo=m-9
    yr=y+1
    mo=m+3
    yr=y
    hour, min, sec =s2hms(secs)
    
    gtime=[yr, mo, d, hour, min, sec]

    return gtime

def julian_day(yyyy,mm,dd,hh):
    """
    Convert date (yyyy, mm, dd, hh) to julian day
    """
    
    IGREG = 15 + 31 * (10+12*1582)
    jy = yyyy
    if (jy < 0):
        jy = jy + 1
    if (mm > 2):
       jm = mm + 1
    else:
       jy = jy - 1
       jm = mm + 13
    
    julian = int(math.floor(365.25*jy)+math.floor(30.6001*jm)+dd+1720995)
    if (dd+31*(mm+12*yyyy) >= IGREG):
       ja = int(0.01*jy)
       julian = julian + 2 - ja + int(0.25*ja)
    
    return julian

def surface_light(time, Lat):
    
    """
    Max light at sea surface
    """
    MAXLIG = 500
    
    # Need day and hour of year
    days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    gtime = gregorian(time)
    
    D = float(gtime[2] + sum(days_in_month[0:int(gtime[1])-1]))
    H = float(gtime[3])
    
    P = math.pi
    TWLIGHT = 5.76
    
    MAXLIG = MAXLIG * 0.1 + MAXLIG * abs(math.sin(math.pi*(D*1.0)/365.0))
    
    """
    Originally the calculation of light as a function of lat, time of day, and time of year was setup for
    sind, cosd of degrees. I converted all to radians using math.radians to suit the python standard.
    """
    DELTA = 0.3979*math.sin(math.radians(0.9856*(D-80)+ 1.9171*(math.sin(math.radians(0.9856*D))-0.98112)))
    H12 = DELTA*math.sin(math.radians(Lat*1.))- math.sqrt(1.-DELTA**2)*math.cos(math.radians(Lat*1.))*math.cos(math.radians(15.0*12.0))
    HEIGHT = DELTA*math.sin(math.radians(Lat*1.))- math.sqrt(1.-DELTA**2)*math.cos(math.radians(Lat*1.))*math.cos(math.radians(15.0*H))
      
    V = math.asin(HEIGHT)
    
    if (V >= 0):                 
        s_light = MAXLIG*(HEIGHT/H12) + TWLIGHT
    elif (V >= -6):
        s_light = ((TWLIGHT - 0.048)/6.)*(6.+V)+.048
    elif (V >= -12):
        s_light = ((0.048 - 1.15e-4)/6.)*(12.+V)+1.15e-4
    elif (V >= -18):
        s_light = (((1.15e-4)-1.15e-5)/6.)*(18.+V)+1.15e-5
    else:
        s_light = 1.15e-5
    print 'func: ', MAXLIG, s_light  
    return s_light
    
    
         
    """
    Test suite: 2440000 should give 1968, 05, 23, 00, 00, 00
    for k in range(24):
        surface_light(2440000+k/24.,42.0)
        print 2440000+k/24.
        
    print julian_day(1968,05,23,00)
    """
    
