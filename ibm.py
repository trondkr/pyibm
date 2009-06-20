import os, sys, string
import netCDF4
from netCDF4 import Dataset
from netCDF4 import num2date
import datetime, types
import numpy as np
import IOwrite

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
    
def init(station):
    """
    init: initializes the reading of files and definition of global variables
    """
    log       = True
    clim      = True
    print 'init',station
    fileNameIn  = station
    
    varlist=['temp','salt','u','v'] #,'nanophytoplankton','diatom','mesozooplankton','microzooplankton','Pzooplankton']
    startDate='15/1/1991'
    endDate='15/7/1991'
    
    """
    Open the netcdf file if it existst.
    """
    cdf = Dataset(fileNameIn)
    """
    Calculate the sigma to meters matrix. This is important as all variables in the netcdf file are stored
    at sigma layers. To be able to convert from sigma to meters we use this function.
    """
    grdSTATION = grd.grdClass(fileNameIn,"STATION")

    """
    Get the time information and find the indices for start and stop data to extract relative to
    the time period wanted. Takes input data:
    startDate="DD/MM/YYYY" and endDate="DD/MM/YYYY" or if none given, finds all date
    """
    getTimeIndices(cdf,grdSTATION,startDate,endDate,clim)
    
    """
    Extract the variables at the given station and store in the grdSTATION object
    """
    IOnetcdf.getStationData(cdf,varlist,grdSTATION,log,clim)
    
    return grdSTATION, clim

def getTimeIndices(cdf,grdSTATION,startDate=None,endDate=None,clim=None):
    """Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """
    ref_date = date.Date()
    ref_date.day=1
    ref_date.month=1
    ref_date.year=1948
    jdref=ref_date.ToJDNumber()
    grdSTATION.jdref=jdref
        
    if clim is False:
        """Get start and stop time stamps from file"""
        startJDFile =int(grdSTATION.time[0])
        endJDFile   =int(grdSTATION.time[-1])
        
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
                    grdSTATION.startIndex=i-1
                    grdSTATION.start_date = date.Date()
                    grdSTATION.start_date =date.DateFromJDNumber(jdstart)
                    FOUND=True
            for i in range(grdSTATION.time.shape[0]):
                if grdSTATION.time[i] +jdref < jdend:
                    FOUND=False
                    continue
                elif grdSTATION.time[i] + jdref >= jdend and FOUND is False:
                    print "Found first time index that fits end point:"
                    print "%s at index %i => %s\n"%(grdSTATION.time[i]+jdref,i,jdend)
                    grdSTATION.endIndex=i+1
                    grdSTATION.end_date = date.Date()
                    grdSTATION.end_date =date.DateFromJDNumber(jdend)
                    
        
                    FOUND=True
        print "Total number of days extracted: %s"%(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])
        grdSTATION.totalDays=(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])
    else:
        print "Using climatological time as clim is set to True"
        d = string.split(startDate,'/')
        start_date = date.Date()
        start_date.day=int(d[0])
        start_date.month=int(d[1])
        start_date.year=int(d[2])
        jdstart=start_date.GetYearDay()
        
        d = string.split(endDate,'/')
        end_date = date.Date()
        end_date.day=int(d[0])
        end_date.month=int(d[1])
        end_date.year=int(d[2])
        jdend=end_date.GetYearDay()
        
        grdSTATION.time=[]
        grdSTATION.time=(cdf.variables["clim_time"][:])*5
        FOUND = False
        
        for i in range(grdSTATION.time.shape[0]):
           
            if grdSTATION.time[i]  < jdstart:
                continue
            elif grdSTATION.time[i]  >= jdstart and FOUND is False:
                print "\nFound first time index that fits start point:"
                print "%s at index %i => %s"%(grdSTATION.time[i],i,jdstart)
                grdSTATION.startIndex=i-1
                grdSTATION.start_date = date.Date()
                grdSTATION.start_date =date.DateFromJDNumber(jdstart)
                FOUND=True
        for i in range(grdSTATION.time.shape[0]):
            if grdSTATION.time[i] < jdend:
                FOUND=False
                continue
            elif grdSTATION.time[i]  >= jdend and FOUND is False:
                print "Found first time index that fits end point:"
                print "%s at index %i => %s\n"%(grdSTATION.time[i],i,jdend)
                grdSTATION.endIndex=i+1
                grdSTATION.end_date = date.Date()
                grdSTATION.end_date =date.DateFromJDNumber(jdend)
                FOUND=True
                
        grdSTATION.jdstart=jdstart
        grdSTATION.jdend=jdend
        print "Total number of days extracted: %s"%(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])
        grdSTATION.totalDays=(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])
        
def ibm(station):
    
    os.system("clear")
    """
    Read the netcdf file and create the input arrays needed by the ibm
    """
    grdSTATION, clim = init(station)
    
    """
    Open output file:
    """
    outputfile=str(grdSTATION.lat[0])+'_'+str(grdSTATION.lon[0])+'_res.nc'
    if os.path.exists(outputfile): os.remove(outputfile)

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
    a = (0.01/3600)*dt;            
    b = -1.3 
    Pe = 0                        
    Ke_larvae = 1
    Ke_predator = 1              
    attenuation_coeff = 0.18
    Ndepths=len(grdSTATION.depth)
    Nhours =24
    Nlarva=50
    Ndays=int(grdSTATION.totalDays)
    Lat = grdSTATION.lat
    
    time=date.Date()
    time.year=grdSTATION.start_date.year
    time.month=grdSTATION.start_date.month
    time.day=grdSTATION.start_date.day
    
    """ julian is the time used to control the larvae, while
    julianFileA and julianFileB are the time stamps of inbetween
    the larvae recides relative to time from file. Temperature is
    interpolated to date julian from julianFileA and julianFileB"""
    
    julian=time.ToJDNumber()
    julianFileA=grdSTATION.time[0]+grdSTATION.jdref
    julianFileB=grdSTATION.time[1]+grdSTATION.jdref
    
    if clim is True:
        julianFileA=grdSTATION.time[0]
        julianFileB=grdSTATION.time[1]
    julianIndex=0
    
    """
    Create array for each individual larvae
    """
    
    larva_array=np.zeros((Nlarva, Nhours, Ndepths, 3),dtype=np.float64) 
    W=np.zeros((Nlarva),dtype=np.float64)
    W_AF=np.zeros((Nlarva),dtype=np.float64)
    S=np.zeros((Nlarva),dtype=np.float64)
    Age=np.zeros((Nlarva),dtype=np.float64)
    L=np.zeros((Nlarva),dtype=np.float64)
    R=np.zeros(13,dtype=np.float64)
    enc=np.zeros(13,dtype=np.float64)
    hand=np.zeros(13,dtype=np.float64)
    pca=np.zeros(13,dtype=np.float64)
    ing=np.zeros(13,dtype=np.float64)
    ingrate=np.zeros((Nlarva),dtype=np.float64)
    Psurvive=np.ones((Nlarva,Nhours),dtype=np.float64)
    Growth=np.zeros((Nlarva, Nhours),dtype=np.float64)
    
    Nprey=1
    grdSTATION.Initialized=False
    grdSTATION.saveIndex=(Ndays,Nlarva,Nprey)
    grdSTATION.Ndays=Ndays
    grdSTATION.Nlarva=Nlarva
    grdSTATION.Nprey=1
    
    prey=0
    
    # Initialize Calanus
    calanus_D  = [0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #(#/ltr)
    calanus_W  = [0.33, 0.49, 1., 1.51,2.09, 2.76, 4.18, 13.24, 23.13, 63.64, 169.58, 276.29, 276.29] #(micrograms)
    calanus_L1 = [0.22, 0.27, 0.4, 0.48, 0.55, 0.61, 0.79, 1.08, 1.38, 1.8, 2.43, 2.11, 2.11] #(length in mm)
    calanus_L2 = [0.1, 0.1, 0.1, 0.15,  0.18, 0.2, 0.22, 0.25, 0.31, 0.41, 0.52, 0.65, 0.65] #(width in mm)

    
    for i in range(Nlarva):
        W[:]    = 0.093 # milligram (5mm larvae)
        W_AF[i] = 0.093
        S[i]    = 0.3*0.06*W[i] # 30% av max mageinnhold
        Age[i] = 0
       
    for day in range(Ndays):
        for ind in range(Nlarva):
            """No behavior, only fixed depth larvae"""
            depth = ind
            depthFOUND=False
            
            for d in range(len(grdSTATION.depth)-1):
                if abs(grdSTATION.depth[d]) < abs(depth) <= abs(grdSTATION.depth[d+1]) and depthFOUND==False: 
                    dz1=1.0-abs((abs(grdSTATION.depth[d])-abs(depth))/(abs(grdSTATION.depth[d])-abs(grdSTATION.depth[d+1])))
                    dz2=1.0-abs((abs(grdSTATION.depth[d+1])-abs(depth))/(abs(grdSTATION.depth[d])-abs(grdSTATION.depth[d+1])))
                    depthIndex1=d; depthIndex2=d+1
                    depthFOUND=True
                    #print 'depth 1:',grdSTATION.depth[d], depth, grdSTATION.depth[d+1],depthIndex1, depthIndex2
                    
                if abs(grdSTATION.depth[0]) > abs(depth) and depthFOUND==False:
                    depthIndex1=d; depthIndex2=d; dz1=0.5; dz2=0.5
                    depthFOUND=True
                    #print 'depth 2:',grdSTATION.depth[d], depth,depthIndex1, depthIndex2
                     
                if abs(grdSTATION.depth[-1]) < abs(depth) and depthFOUND==False:
                    depthIndex1=-1; depthIndex2=-1; dz1=0.5; dz2=0.5
                    depthFOUND=True
                    #print 'depth 3:',grdSTATION.depth[-1], depth,depthIndex1, depthIndex2
                    
            
            # TODO: Fix so that you can use any timestep. Do this by making the loop
            # count on Nsteps instead of Nhours, where Nsteps*dt=Nhours.
            
            aveLight=0.0
            for hour in range(Nhours):
                L[ind] = np.exp(2.296 + 0.277*np.log(W[ind]) - 0.005128*np.log(W[ind])**2)
                
                currentDate=IOtime.julian_to_date(julian, hour)
                julian     =IOtime.julian_day(currentDate[0],currentDate[1],currentDate[2],hour)
                
                """Calculate weights to use on input data from file"""
                dwA = abs(julian) - abs(julianFileA)
                dwB = abs(julianFileB) - abs(julian)
                
                s_light = IOlight.surface_light(julian,grdSTATION.lat,hour);
               
                Eb = s_light*np.exp(attenuation_coeff*(-depth));
                aveLight=aveLight + Eb
                """
                Mouthsize of larvae. The larvae can only capture calanus of sizes less
                than the larval mouthsize. Folkvord et al.
                """
                m = np.exp (-3.27+1.818*np.log(L[ind])-0.1219*(np.log(L[ind]))**2.)
                
                meta = dt*2.38e-7*np.exp(0.088*7)*((W[ind]*1000)**(0.9)*0.001)
                if Eb > 0.001:
                    if L[ind] > 5.5:
                        meta = (2.5*meta)
                    else:
                        meta = (1.4*meta);
                """Interpolate the values of temp, salt, u and v velocity in time to current julian date"""
                Tdata=((grdSTATION.data[julianIndex,depthIndex1,0])*
                    (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,0])*
                    (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,0])*
                    (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,0])*
                    (dwB/(dwA+dwB)))*dz2
                Sdata=((grdSTATION.data[julianIndex,depthIndex1,1])*
                    (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,1])*
                    (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,1])*
                    (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,1])*
                    (dwB/(dwA+dwB)))*dz2
                Udata=((grdSTATION.data[julianIndex,depthIndex1,2])*
                    (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,2])*
                    (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,2])*
                    (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,2])*
                    (dwB/(dwA+dwB)))*dz2
                Vdata=((grdSTATION.data[julianIndex,depthIndex1,3])*
                    (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,3])*
                    (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,3])*
                    (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,3])*
                    (dwB/(dwA+dwB)))*dz2
                
                assi = 0.8*(1.0 - 0.4*np.exp(-0.002*(W[ind]/mm2m-30.0)))
                GR = sec2day*np.log(0.01*(1.08 + 1.79*Tdata - 0.074*Tdata*np.log(W[ind])
                                          - 0.0965*Tdata*np.log(W[ind])**2 + 0.0112*Tdata*np.log(W[ind])**3) + 1)*dt  
                GR_gram = (np.exp(GR) - 1)*W[ind]
                               
               
                for j in range(13):
                    Ap_calanus = 0.75*calanus_L1[j]*mm2m*calanus_L2[j]*mm2m
                    R[j] = IOlight.get_perception_distance(attenuation_coeff,Ke_larvae,Ap_calanus,Eb)
                  

                for j in range(13):
                    hand[j] = 0.264*10**(7.0151*(calanus_L1[j]/L[ind])) # Walton 1992
                    enc[j] = (0.667*np.pi*(R[j]**3.)*f + np.pi*(R[j]**2.)*np.sqrt(calanus_L1[j]**2.+ 2.*omega**2.)*f*tau) * calanus_D[j]* ltr2mm3
                    pca[j] = enc[j]*max(0.0,min(1.0,-16.7*(calanus_L1[j]/L[ind]) + 3.0/2.0))
                    ing[j] = dt*pca[j]*calanus_W[j]*micro2m / (1 + hand[j]);
                  
                W_old = W[ind]
                S_old = S[ind]
                S[ind] = min(gut_size*W[ind],S_old + sum(ing[:]))
                ingrate[ind] = (S[ind] - S_old)/W[ind]
                W[ind] = W[ind] + min(GR_gram + meta,S[ind]*assi) - meta #- 0.1*act*abs(Init(l,6) - Init(l,6))/(L(l)*dt)*meta;
                S[ind] = S[ind] - ((W[ind] - W_old) + meta)/assi # + 0.1*act*abs(Init(l,6) - Init(l,6))/(L(l)*dt)*meta)/assi;
                
                W_AF[ind] = W_AF[ind] + (np.exp(GR)-1)*W_AF[ind]
         
                L[ind] = min(L[ind],np.exp(2.296 + 0.277*np.log(W[ind]) - 0.005128*np.log(W[ind])**2))
                Ap_larvae = 0.1*L[ind]**2*mm2m**2
                Rpisci = IOlight.get_perception_distance(attenuation_coeff,Ke_predator,Ap_larvae,Eb)
              
                Psurvive[ind,hour] = np.exp(-a*(L[ind]**b) - C2*(1-Pe)*Rpisci**2)*float(Psurvive[ind,hour])
                """Sum up the total SGR over 24 hours"""
                Growth[ind, hour]=((W[ind]-W_old)/W_old)*100.0
                
                
            if depth <= Nlarva-1:
                if currentDate[1]<10:
                    print_mm='0'+str(currentDate[1])
                else:
                    print_mm=str(currentDate[1])
                if time.day<10:
                    print_dd='0'+str(currentDate[2])
                else:
                    print_dd=str(currentDate[2])
                if clim is False:   
                    print_date=str(currentDate[0])+'-'+print_mm+'-'+str(print_dd)+'T%2i:00'%(hour)
                else:
                    print_date='CLIM-'+print_mm+'-'+str(print_dd)+'T%2i:00'%(hour)
            
            """Calculate 24 hour SGR """        
            max_g = 1.08 + 1.79 * Tdata - 0.074 * Tdata*np.log(W[ind]) - 0.0965 * Tdata *((np.log(W[ind]))**2) + 0.0112 * Tdata *((np.log(W[ind])**3.))
          
          
            """Store results in grdSTATION object and use that object to write to netCDF4 file"""                  
            grdSTATION.larvaTime=julian
            grdSTATION.larvaDepth=-depth
            grdSTATION.larvaWgt=W[ind]
            grdSTATION.larvaSgr=sum(Growth[ind])
            grdSTATION.larvaSgrAF=max_g
                
            grdSTATION.larvaAF=W_AF[ind]
            grdSTATION.larvaPsur=np.exp(-sum(Psurvive[ind,:]))*100.
            grdSTATION.larvaTdata=Tdata
            grdSTATION.larvaAveLight=aveLight/24.0
            
            IOwrite.writeStationFile(grdSTATION,day,depth,prey,outputfile)
                
            if ind==(Nlarva-1): 
                julian=julian+1
                print '%3.1f\t %s\t %4.2f'%(grdSTATION.lat, print_date,julian)
            if julian >=julianFileB and julian < grdSTATION.jdend:
               # print 'reset time julian',julianFileA,julian,julianFileB  ,  grdSTATION.jdend
                julianIndex=julianIndex+1
                if clim is False:
                    julianFileA=grdSTATION.time[julianIndex]+grdSTATION.jdref
                    julianFileB=grdSTATION.time[julianIndex+1]+grdSTATION.jdref
                else:
                    julianFileA=grdSTATION.time[julianIndex]
                    julianFileB=grdSTATION.time[julianIndex+1]
            if julian == grdSTATION.jdend:
                #print 'Finished running'
                continue
            #  Resetting values for larvae.
           
            if ind==(Nlarva-1): 
                for i in range(Nlarva):
                    W[:]    = 0.093 # milligram (5mm larvae)
                    W_AF[i] = 0.093
                    S[i]    = 0.3*0.06*W[i] # 30% av max mageinnhold
                    Age[i] = 0
                    Psurvive[i,:]=0.0
            
                                                             
            hour=0
        
if __name__=="__main__":
    import psyco
    try:
        import psyco
        psyco.log()
        psyco.full(memory=100)
        psyco.profile(0.05, memory=100)
        psyco.profile(0.2)
    except ImportError:
        pass
    
    stations=["/Users/trond/Projects/arcwarm/SODA/soda2roms/station_lat_41.5423_lon_-132.9234.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/station_lat_43.4111_lon_-50.4321.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/station_lat_44.5001_lon_-54.3801.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/station_lat_63.3501_lon_1.5301.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/station_lat_64.7201_lon_21.5101.nc"]
    #stations=["/Users/trond/Projects/arcwarm/SODA/soda2roms/station_lat_63.3501_lon_1.5301.nc"]

    
    for station in stations:
        
        ibm(station)
  