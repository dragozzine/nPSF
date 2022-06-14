import numpy as np
import pandas as pd
from astroquery.jplhorizons import Horizons
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
from astropy.coordinates import GeocentricMeanEcliptic
from astropy.coordinates import HeliocentricMeanEcliptic
from astropy.coordinates import HeliocentricTrueEcliptic
from astropy.coordinates import GeocentricTrueEcliptic
from astropy.coordinates import BarycentricTrueEcliptic
from astropy.coordinates import BarycentricMeanEcliptic
from astropy import coordinates 
import astropy
import sys
import commentjson as json
import gc


'''
    NAME:
         convert_to_primary_centric
         
    PURPOSE:
         This function takes a parameter Dataframe in RA/DEC, and converts it to Latitude 
         and Longitude, while also converting the dates to Primary-Centric Julian dates
         
    CALLING SEQUENCE:
         convert_to_primary_centric(paramsDF, objectName)
   
    INPUTS
          paramsDF - A dataframe of the observed positional data of the KBO in question
          objectName - The name of the object being observed (needed for the Horizons function)
   
    OUTPUTS:
          None. Just makes plots currently.
'''


#From latlon_transform.py
def convert_to_primary_centric(paramsDF, objectNames, numobjects, resultspath, sample_num):
     #Current column names are just descriptive, not representative of final product column names
    updatedDF = pd.DataFrame(columns = ['time'])
    forsigsDF = pd.DataFrame(columns = ['NaN_list'])
    
    #Convert the dates into a Julian date format
    date = paramsDF['time'].dropna()
    dateList = []
    for i in date:
            jd = Time(i,format='jd')
            dateList.append(jd)
    #print(dateList)
    # Make new goecentric position file
    ##import mm_make_geo_pos
    ##geo = mm_make_geo_pos.mm_make_geo_pos(objectNames[0], dateList, runprops = None, synthetic = True)
    ##geo.to_csv("geocentric_" + objectNames[0] + "_position.csv")
        
    #print(dateList)
    #Get the Horizons data for the object at the times it was observed
    primary = Horizons(id=objectNames[0],location="geocentric",epochs=dateList)
    #primary = Horizons(id=objectNames[0],location=None,epochs=dateList)
    
    updatedDF['time'] = paramsDF['time']-primary.vectors(aberrations = 'astrometric')['lighttime']

    #Pull all data from csv file
    #RA_Prim = np.array(paramsDF['RA-Primary'])
    #DEC_Prim = np.array(paramsDF['DEC-Primary'])
    RA_Prim = np.array(primary.ephemerides()['RA'][:])
    DEC_Prim = np.array(primary.ephemerides()['DEC'][:])
    print(RA_Prim, DEC_Prim)
    
    num = 1
    for i in range(len(objectNames)-1):
    
        deltaRA_1 = np.array(paramsDF['Delta-RA_'+objectNames[i+1]]).astype(np.float)
        deltaDEC_1 = np.array(paramsDF['Delta-DEC_'+objectNames[i+1]]).astype(np.float)
    
        RA_1_err = np.array(paramsDF['Delta-RA_'+objectNames[i+1]+'-err']).astype(np.float)
        DEC_1_err = np.array(paramsDF['Delta-DEC_'+objectNames[i+1]+'-err']).astype(np.float)

        RA_1 = RA_Prim+deltaRA_1/3600/np.cos(DEC_Prim*u.degree)
        DEC_1 = DEC_Prim + deltaDEC_1/3600
    
        ra_err = np.zeros((len(RA_1), int(sample_num)))
        dec_err = np.zeros((len(DEC_1), int(sample_num)))
        
        #Here we create the randomnly distributed ra and dec errors
        #print(RA_1, len(RA_1), RA_1_err, len(RA_1_err))
        print(len(RA_1),len(RA_1_err))
        
        #for k in tqdm(range(len(RA_1))):
            #plt.figure(k)
            #for j in range(sample_num):
                #ra_err[k][j] = np.random.normal(RA_1[k]*3600, abs(RA_1_err[k]))/3600
                #dec_err[k][j] = np.random.normal(DEC_1[k]*3600, abs(DEC_1_err[k]))/3600
            #plt.scatter(ra_err[k],dec_err[k],s=10)
            
        for k in tqdm(range(len(RA_1))):
            ra_err[k] = np.random.normal(RA_1[k]*3600, abs(RA_1_err[k]), (int(sample_num)))/3600
            dec_err[k] = np.random.normal(DEC_1[k]*3600, abs(DEC_1_err[k]), (int(sample_num)))/3600
        
    #Essentially we define where the object is in our RA/DEC coordinate system. ICRS is the system our coordinates are in.
        dist = primary.vectors(aberrations = 'astrometric')['range']

        firstC = SkyCoord(ra=RA_1*u.degree, dec=DEC_1*u.degree, frame='gcrs', obstime = dateList, distance = dist)
        primC = SkyCoord(ra=RA_Prim*u.degree, dec=DEC_Prim*u.degree, frame='gcrs', obstime = dateList, distance = dist)
        firstEcl = firstC.transform_to(GeocentricTrueEcliptic(equinox='J2000'))
        primEcl = primC.transform_to(GeocentricTrueEcliptic(equinox='J2000'))
    
        Lat_Prim = primEcl.lat.degree
        Long_Prim = primEcl.lon.degree
        #print(Lat_Prim, Long_Prim)
    
        Lat_1 = firstEcl.lat.degree
        Long_1 = firstEcl.lon.degree
        #print(Lat_1, Long_1)
                
        Lat_err = np.zeros(len(ra_err))
        Long_err = np.zeros(len(dec_err))
        
        #transform all of the randomnly distributed errors
        print("Begin error transform")
        coord_sky = SkyCoord(ra=ra_err*u.degree, dec=dec_err*u.degree, frame='gcrs', obstime = dateList[0], distance = dist[0]*u.AU,unit=(u.deg,u.deg))
        print("Begin GeocentricTrueEcliptic transform")
        transformed_coord = coord_sky.transform_to(GeocentricTrueEcliptic(equinox='J2000'))
        print("Begin lat degree transform")
        Lat_err_arr = transformed_coord.lat.degree
        print("Begin lon degree transform")
        Long_err_arr = transformed_coord.lon.degree

        for j in range(len(Lat_err_arr)):
            Lat_err[j] = np.std(Lat_err_arr[j])
            Long_err[j] = np.std(Long_err_arr[j])
        
        # clear astroquery memory, necessary for running more than 2psfs.
        del coord_sky
        del transformed_coord
        del Lat_err_arr
        del Long_err_arr
        gc.collect()
        
        #transform all of the randomnly distributed errors
        #for j in range(len(ra_err)):
            #Lat_err_arr = np.zeros(len(ra_err[0]))
            #Long_err_arr = np.zeros(len(dec_err[0]))
            #for k in range(len(ra_err[0])):
                #coord_sky = SkyCoord(ra=ra_err[j][k]*u.degree, dec=dec_err[j][k]*u.degree, frame='gcrs', obstime = dateList[0], distance = dist[0]*u.AU,unit=(u.deg,u.deg))
                #transformed_coord = coord_sky.transform_to(GeocentricTrueEcliptic(equinox='J2000'))
                #Lat_err_arr[k] = transformed_coord.lat.degree
                #Long_err_arr[k] = transformed_coord.lon.degree
            #Lat_err[j] = np.std(Lat_err_arr)
            #Long_err[j] = np.std(Long_err_arr)
            
        #print(np.array_equal(Laterr, Lat_err))
        #print(np.array_equal(Longerr,Long_err))
        
        DeltaLat_1 = (Lat_1-Lat_Prim)*3600
        DeltaLong_1 = (Long_1-Long_Prim)*np.cos(Lat_Prim*u.degree)*3600
        
        # Ask Dallin about this code!!!
        #DeltaLong_1 = (Long_1-Long_Prim)*np.cos(Lat_Prim*u.degree)*3600
        #print(DeltaLong_1)
        #DeltaLong_1 = (Long_1-Long_Prim)*np.cos(np.deg2rad(Lat_Prim))*3600 
        #print(DeltaLong_1)
    
        Lat_1_err_arc = (Lat_err)*3600
        Long_1_err_arc = (Long_err)*3600
    
    
        #if i == 0:
            #updatedDF['Lat_Prim'] = Lat_Prim
            #updatedDF['Long_Prim'] = Long_Prim
    
        updatedDF['DeltaLat_'+objectNames[i+1]] = DeltaLat_1
        updatedDF['DeltaLong_'+objectNames[i+1]] = DeltaLong_1

        updatedDF['DeltaLat_'+objectNames[i+1]+'_err'] = Lat_1_err_arc
        updatedDF['DeltaLong_'+objectNames[i+1]+'_err'] = Long_1_err_arc
        
        forsigsDF['lat'+ str(numobjects - num)] = Lat_1
        forsigsDF['long'+ str(numobjects - num)] = Long_1
        forsigsDF['dlat'+ str(numobjects - num)] = DeltaLat_1
        forsigsDF['dlong'+ str(numobjects - num)] = DeltaLong_1
        num -= 1        

        
    #create a single line output obs_df file for multimoon
    obsDF = pd.DataFrame(columns = list(updatedDF.columns), index = ["0"])
    obsDF['time'] = updatedDF['time'].iloc[0]
    
    for i in range(len(updatedDF.columns)-1):
        num = updatedDF.iloc[:,i+1]

        median = np.percentile(num,50, axis = None)
        
        obsDF.iloc[:,i+1] = median

    obsDF.insert(1, 'Lat_Prim', Lat_Prim) 
    obsDF.insert(2, 'Long_Prim', Long_Prim)
    print(obsDF)

    filename = objectNames[0]+'_obs_df.csv'
    obsDF.to_csv(resultspath + "/" + filename)
    
    return forsigsDF
    

   