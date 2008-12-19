import math, datetime

"""
Functions for conversion between gregorian calendar, julian day, and regular date format.

Converted from Fortran to Python on 10.06.2008 by Trond Kristiansen
Trond.Kristiansen@imr.no

Functions in module:

1. [hour, min, sec]=s2hms(secs)
2. [yyyy, month, day, hour, min, sec]=gregorian(julian)
3. [yyyy, month, day, hour, 00, 00]=julian_to_date(julian)
4. [julian]=julian_day(yyyy,mm,dd,hh)
"""

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 10)
__modified__ = datetime.datetime(2008, 6, 10)
__version__  = "1.1"
__status__   = "Production"

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
    print yr, mo, d, hour, sec
    gtime=[yr, mo, d, hour, min, sec]

    return gtime

def julian_to_date(julian, hour):
    IGREG = 2299161
    if (julian >= IGREG):
       x = ((julian-1867216)-0.25) / 36524.25
       ja = julian + 1 + int (x) - int (0.25*x)
    else:
       ja = julian
    
    jb = ja + 1524
    jc = int(6680+((jb-2439870)-122.1)/365.25)
    jd = int(365*jc+(0.25*jc))
    je = int((jb-jd)/30.6001)
    
    dd = jb - jd - int(30.6001*je)
    mm = je - 1
    if (mm > 12):
        mm = mm - 12
    yyyy = jc - 4715
    if (mm > 2):
        yyyy = yyyy - 1
    if (yyyy <= 0):
        yyyy = yyyy - 1
    
    day = dd
    year = yyyy
    month = mm
    gtime=[yyyy, month, day, hour, 00, 00]
    
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