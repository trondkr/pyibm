import os, sys, string
import netCDF4
from netCDF4 import Dataset
from netCDF4 import num2date
import datetime, types, math
import numpy as np

"""
Import modules created as part of project
"""

""" Get self made modules"""
dir='/Users/trond/Projects/arcwarm/SODA/soda2roms'

if os.path.isdir(dir):
    sys.path.append(dir)


import grd
import IOverticalGrid
import IOtime
import IOlight
import date
import IOnetcdf

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 9)
__modified__ = datetime.datetime(2008, 6, 9)
__version__  = "1.1.5"
__status__   = "Production"

def help():
    """
    Print the helpful module docstring.
    help() -> None
    """
    
    print __doc__
    return
    
def init():
    """
    init: initializes the reading of files and definition of global variables
    """
    log       = True
    fileNameIn  = "/Users/trond/Projects/arcwarm/SODA/soda2roms/station_lat_41.5423_lon_-66.5233.xyz"
    
    varlist=['temp','salt','u','v'] #,'nanophytoplankton','diatom','mesozooplankton','microzooplankton','Pzooplankton']
    startDate='20/2/1991'
    endDate='10/10/1994'
    
    """
    Open the netcdf file if it existst.
    """
    cdf = Dataset(fileNameIn)
    """
    Calculate the sigma to meters matrix. This is important as all variables in the netcdf file are stored
    at sigma layers. To be able to convert from sigma to meters we use this function.
    """
    grdSTATION = grd.grdClass(fileNameIn,"STATION")
    #IOverticalGrid.get_z_levels(grdSTATION)

    """
    Get the time information and find the indices for start and stop data to extract relative to
    the time period wanted. Takes input data:
    startDate="DD/MM/YYYY" and endDate="DD/MM/YYYY" or if none given, finds all date
    """
    getTimeIndices(cdf,grdSTATION,startDate,endDate)
    
    """
    Extract the variables at the given station and store in the grdSTATION object
    """
    IOnetcdf.getStationData(cdf,varlist,grdSTATION,log)
    
    return grdSTATION

def getTimeIndices(cdf,grdSTATION,startDate=None,endDate=None):
    """Get start and stop time stamps from file"""
    startJDFile =int(grdSTATION.time[0])
    endJDFile   =int(grdSTATION.time[-1])
    """
    Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """
    ref_date = date.Date()
    ref_date.day=1
    ref_date.month=1
    ref_date.year=1948
    jdref=ref_date.ToJDNumber()
    
    """Figure out the date from the JD number. I have consistently used JD numbers
    relative to 01/01/1948 and you therefore have to add that jdref number
    to calculate the correct date using the date module."""
 
    startDateFile=date.DateFromJDNumber(startJDFile+jdref)
    endDateFile=date.DateFromJDNumber(endJDFile+jdref)
   
    print "\nStation contains data for the time period:"
    print "%s to %s"%(startDateFile, endDateFile)
    
    """Now find the Juian dates of the time period you have asked for and see
    if it exists in the file. If so, the find the indices that correspond to the
    time period.
    """

    d = string.split(startDate,'/')
    start_date = date.Date()
    start_date.day=int(d[0])
    start_date.month=int(d[1])
    start_date.year=int(d[2])
    jdstart=start_date.ToJDNumber()

    d = string.split(endDate,'/')
    end_date = date.Date()
    end_date.day=int(d[0])
    end_date.month=int(d[1])
    end_date.year=int(d[2])
    jdend=end_date.ToJDNumber()
    
    if jdstart < startJDFile+jdref:
        print "Start time preceeds the earliest time stamp found in file"
        print "Required in File: %i/%i/%i"%(int(end_date.day),int(end_date.month),int(end_date.year))
        print "Actually in File: %i/%i/%i"%(int(endDateFile.day),int(endDateFile.month),int(endDateFile.year))
        exit()
        
    if jdend > endJDFile+jdref:
        print "End time proceeds the latest time stamp found in file"
        print "Required in File: %i/%i/%i"%(int(start_date.day),int(start_date.month),int(start_date.year))
        print "Actually in File: %i/%i/%i"%(int(startDateFile.day),int(startDateFile.month),int(startDateFile.year))
        exit()
        
    if jdend < endJDFile+jdref and jdstart > startJDFile+jdref:
        print "Time period to extract was found within the time period available in the file..."
        print "--> %s - %s"%(startDate,endDate)
        
        for i in range(grdSTATION.time.shape[0]):
            if grdSTATION.time[i] +jdref < jdstart:
                FOUND=False
                continue
            elif grdSTATION.time[i] + jdref >= jdstart and FOUND is False:
                print "\nFound first time index that fits start point:"
                print "%s at index %i => %s"%(grdSTATION.time[i]+jdref,i,jdstart)
                grdSTATION.startIndex=i
                FOUND=True
        for i in range(grdSTATION.time.shape[0]):
            if grdSTATION.time[i] +jdref < jdend:
                FOUND=False
                continue
            elif grdSTATION.time[i] + jdref >= jdend and FOUND is False:
                print "Found first time index that fits end point:"
                print "%s at index %i => %s\n"%(grdSTATION.time[i]+jdref,i,jdend)
                grdSTATION.endIndex=i
                FOUND=True
    print "Total number of days extracted: %s"%(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])
    grdSTATION.totalDays=(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])
    
def ibm():
    os.system("clear")
    """
    Read the netcdf file and create the input arrays needed by the ibm
    """
    grdSTATION = init()
    
    """
    Define global variables
    """
    dt = 3600;                   
    sec2day = 1.0/86400.0
    gut_size = 0.06
    stomach_threshold = 0.3
    tau = 2.0                     
    omega = 0.0
    f = 0.43                     
    pi = 3.14159265
    mm2m = 0.001
    ltr2mm3 = 1e-6
    micro2m = 0.001
    C2 = 0.05                     
    act = 1
    a = 0.01/3600*dt;            
    b = -1.3 
    Pe = 0                        
    Ke_larvae = 1
    Ke_predator = 1              
    attenuation_coeff = 0.18
    Ndepths=len(grdSTATION.depth)
    Nhours =24
    Nlarva=29
    Ndays=grdSTATION.totalDays
    Lat = 41.0
    
    # TODO: Init year, month, day should be from netcdf file
    yyyy=2007
    mm=1
    dd=1
    julian = IOtime.julian_day(yyyy,mm,dd,0)
    """
    Create array for each individual larvae
    """
    
    larva_array=np.zeros((Nlarva, Nhours, Ndepths, 3),float) 
    W=np.zeros((Nlarva),float)
    W_AF=np.zeros((Nlarva),float)
    S=np.zeros((Nlarva),float)
    Age=np.zeros((Nlarva),float)
    L=np.zeros((Nlarva),float)
    R=np.zeros(13,float)
    enc=np.zeros(13,float)
    hand=np.zeros(13,float)
    pca=np.zeros(13,float)
    ing=np.zeros(13,float)
    ingrate=np.zeros((Nlarva),float)
    Psurvive=np.ones((Nlarva),float)
    Growth=np.zeros((Nlarva, Nhours),float)
    
    
    # Initialize Calanus
    calanus_D  = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #(#/ltr)
    calanus_W  = [0.33, 0.49, 1., 1.51,2.09, 2.76, 4.18, 13.24, 23.13, 63.64, 169.58, 276.29, 276.29] #(micrograms)
    calanus_L1 = [0.22, 0.27, 0.4, 0.48, 0.55, 0.61, 0.79, 1.08, 1.38, 1.8, 2.43, 2.11, 2.11] #(length in mm)
    calanus_L2 = [0.1, 0.1, 0.1, 0.15,  0.18, 0.2, 0.22, 0.25, 0.31, 0.41, 0.52, 0.65, 0.65] #(width in mm)

    for i in range(Nlarva):
        W[:]    = 0.093 # milligram (5mm larvae)
        W_AF[i] = 0.093
        S[i]    = 0.3*0.06*W[i] # 30% av max mageinnhold
        Age[i] = 0
       
    for day in range(Ndays):
        julian=julian+1
        for ind in range(Nlarva):
                depth = abs(z_r[day,ind,eta,xi])
                
                for hour in range(Nhours):
                    L[ind] = math.exp(2.296 + 0.277*math.log(W[ind]) - 0.005128*math.log(W[ind])**2)
                    
                    gtime=IOtime.julian_to_date(julian, hour)
                    
                    dd=gtime[2]
                    mm=gtime[1]
                    yyy=gtime[0]
                    julian = IOtime.julian_day(yyyy,mm,dd,hour)
                   
                    s_light = IOlight.surface_light(2440000, Lat, hour);
                   
                    Eb = s_light*math.exp(attenuation_coeff*(-depth));
                   
                    """
                    Mouthsize of larvae. The larvae can only capture calanus of sizes less
                    than the larval mouthsize. Folkvord et al.
                    """
                    m = math.exp (-3.27+1.818*math.log(L[ind])-0.1219*(math.log(L[ind]))**2.)
                    
                    meta = dt*2.38e-7*math.exp(0.088*7)*((W[ind]*1000)**(0.9)*0.001)
                    if Eb > 0.001:
                        if L[ind] > 5.5:
                            meta = (2.5*meta)
                        else:
                            meta = (1.4*meta);
                    
                    T=var_array[day,ind,0]
                    
                    assi = 0.8*(1.0 - 0.4*math.exp(-0.002*(W[ind]/mm2m-30.0)))
                    GR = sec2day*math.log(0.01*(1.08 + 1.79*T - 0.074*T*math.log(W[ind]) - 0.0965*T*math.log(W[ind])**2 + 0.0112*T*math.log(W[ind])**3) + 1)*dt  
                    GR_gram = (math.exp(GR) - 1)*W[ind]
                                   
                   
                    for j in range(13):
                        Ap_calanus = 0.75*calanus_L1[j]*mm2m*calanus_L2[j]*mm2m
                        R[j] = IOlight.get_perception_distance(attenuation_coeff,Ke_larvae,Ap_calanus,Eb)
                      

                    for j in range(13):
                        hand[j] = 0.264*10**(7.0151*(calanus_L1[j]/L[ind])) # Walton 1992
                        enc[j] = (0.667*math.pi*(R[j]**3.)*f + math.pi*(R[j]**2.)*math.sqrt(calanus_L1[j]**2.+ 2.*omega**2.)*f*tau) * calanus_D[j]* ltr2mm3*(var_array[day,ind,3]*100.0)
                        pca[j] = enc[j]*max(0.0,min(1.0,-16.7*(calanus_L1[j]/L[ind]) + 3.0/2.0))
                        ing[j] = 3600.*pca[j]*calanus_W[j]*micro2m / (1 + hand[j]);
                      
                    W_old = W[ind];
                    S_old = S[ind];
                    S[ind] = min(gut_size*W[ind],S_old + sum(ing[:]))
                    ingrate[ind] = (S[ind] - S_old)/W[ind]
                    W[ind] = W[ind] + min(GR_gram + meta,S[ind]*assi) - meta #- 0.1*act*abs(Init(l,6) - Init(l,6))/(L(l)*dt)*meta;
                    S[ind] = S[ind] - ((W[ind] - W_old) + meta)/assi # + 0.1*act*abs(Init(l,6) - Init(l,6))/(L(l)*dt)*meta)/assi;
                    
                    W_AF[ind] = W_AF[ind] + (math.exp(GR)-1)*W_AF[ind]
             
                    L[ind] = min(L[ind],math.exp(2.296 + 0.277*math.log(W[ind]) - 0.005128*math.log(W[ind])**2))
                    Ap_larvae = 0.1*L[ind]**2*mm2m**2
                    Rpisci = IOlight.get_perception_distance(attenuation_coeff,Ke_predator,Ap_larvae,Eb)
                  
                    Psurvive[ind] = math.exp(-a*(L[ind]**b) - C2*(1-Pe)*Rpisci**2)*float(Psurvive[ind])
                    
                    
                    Growth[ind, hour]=((W[ind]-W_old)/W_old)*100.0
                if depth <= 60:
                    if mm<10:
                        print_mm='0'+str(mm)
                    else:
                        print_mm=str(mm)
                    if dd<10:
                        print_dd='0'+str(dd)
                    else:
                        print_dd=str(dd)
                    print_date=str(yyyy)+'-'+print_mm+'-'+str(print_dd)+'T00:00'
                    print '%s\t %f\t %f\t %s'%(print_date, -depth, sum(Growth[ind]), var_array[day,ind,6])
                #  Resetting values for larvae.
                for i in range(Nlarva):
                    W[:]    = 0.093 # milligram (5mm larvae)
                    W_AF[i] = 0.093
                    S[i]    = 0.3*0.06*W[i] # 30% av max mageinnhold
                    Age[i] = 0
                #print '\n------day %i ------------------------'%(day)
                #print 'Survival prob: %10.10f '%(Psurvive[ind])
                #print 'Larvae at depth : %s'%(depth)
                #print 'Weight %f and length %f of individual %i '%(W[ind], L[ind], ind)
                #print 'Growth rate %f'%(sum(Growth[ind,:]))
                                                                 
                     
if __name__=="__main__":
    ibm()
  