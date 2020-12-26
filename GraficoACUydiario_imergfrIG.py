#################################################################################
#script para graficar los datos de precipitacion estimada por IMERG
#Adaptado por Martin Iglesias Github SudestadaARG 

#IMERG Integrated Multi-satellitE Retrievals for GPM
#https://arset.gsfc.nasa.gov/sites/default/files/water/Brazil_2017/Day4/S7P1-span.pdf
#https://github.com/SaulMontoya/Representacion-Interactiva-de-Precipitacion-IMERG-con-Python---Caso-Colombia
########################################################################
import h5py
import numpy as np
import os
from datetime import datetime  
from datetime import timedelta
import argparse


# En caso de querer un elemento por linea de comando:
#parser = argparse.ArgumentParser(description='Year Month Day File')
#parser.add_argument('Year',type=int)
#parser.add_argument('Month',type=int)
#parser.add_argument('Day',type=int)
#parser.add_argument('File',type=str)

#YY= args.Year
#MM=args.Month
#DD= args.Day
#Datafile= arg.File
Fecha='20151219'
Datafile='3B-HHR.MS.MRG.3IMERG.20151219-S000000-E002959.0000.V06B.HDF5' #tomo el primer dia del periodo 
f = h5py.File('/data/miglesias/IMERGFR/'+Fecha+'/'+Datafile, 'r')
a_group_key = list(f.keys())[0]
data = f[a_group_key]	# con esto te quedas solo con "Grid", y adentro tenes las siguientes variables:
data.keys()
# [u'nv', u'lonv', u'latv', u'time', u'lon', u'lat', u'time_bnds', u'lon_bnds', u'lat_bnds', u'precipitationCal', u'precipitationUncal', u'randomError', u'HQprecipitation', u'HQprecipSource', u'HQobservationTime', u'IRprecipitation', u'IRkalmanFilterWeight', u'probabilityLiquidPrecipitation', u'precipitationQualityIndex']


# Lo que nos interesa es precipitationCal
#https://pmm.nasa.gov/sites/default/files/document_files/IMERG_doc.pdf
#The *calibrated precipitation field*, data field precipitationCal, provides global Level 3
#(0.1grx0.1gr-gridded) half-hourly listing in 3IMERGHH data files of the calibrated data described in final 'post-processing'. In 3IMERGM data files it provides the global Level 3 (0.1grx0.1gr- gridded) monthly SG combination data described in 'SG combination'. In the Early and Late Runs, the field provides data with month/location-varying climatological calibration to #the Final Run.
#estimates from all constellation members in the microwave-only  precipitation field (HQprecipitation) and the complete precipitation fields (precipitationCal, precipitationUncal) over #the fully global domain (90grN-S). 
#!Note: microwave estimates over snowy-icy surface types are not masked out in HQprecipitation, but are for precipitationCal and precipitationUncal (as was done in V05). 




# Saco las variables que me interesan:
lat = np.array(data.get('lat'))[:]	# dim 1800 -> resolucion 0.1 grados
lon = np.array(data.get('lon'))[:]	# dim 3600 -> resolucion 0.1 grados
pp = np.array(data.get('precipitationCal'))
pp[pp < 0] =0

lat=lat[499:651]	# de -40.05 a -24.95
lon=lon[1149:1301]	# de -65.05 a -49.95
lon, lat = np.meshgrid(lon, lat)


#########################################################################
#Acumulado Periodo (se hace una vez)
# para abrir varios_:
DataPath = '/data/miglesias/IMERGFR/Acu_2015_12_19-23/'
Filelist = sorted(os.listdir(DataPath))
ppest=np.ma.zeros([len(Filelist),152,152])	# [tiempos,lon,lat]
for nfile in np.arange(0,len(Filelist),1):
	f = h5py.File(DataPath + Filelist[nfile], 'r')
	a_group_key = list(f.keys())[0]
	data = f[a_group_key]
	pp = np.array(data.get('precipitationCal'))
	ppest[nfile,:,:] = pp[0,1149:1301,499:651] #aca estan en lon lat todavia
	ppest[ppest < 0] =0

ppacum = np.sum(ppest, axis=0)	# sumo todos los tiempos, acumulado total del periodo
ppacum = ppacum*0.5   # multiplico por 0.5 para obtener mm/h
ppacum = np.transpose(ppacum)	# necesito [lat, lon] y los datos vienen en [lon, lat]
# en el script de NCL de Cyn tambien intercambia el orden de lat, lon mediante pp = data(lat|:,lon|:)

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm

plt.figure(figsize=(10,10))
intervalos = [1,5,10,25,50,100,150,200,250]
latcorners = ([-40,-25])
loncorners = ([-65,-50])
m = Basemap(projection='cyl',llcrnrlat=latcorners[0],urcrnrlat=latcorners[1],llcrnrlon=loncorners[0],urcrnrlon=loncorners[1],resolution='l')
plt.contourf(lon,lat,ppacum,intervalos,cmap=plt.cm.viridis_r,extend="both")
cbar = plt.colorbar(ticks=intervalos,shrink=0.5)
cbar.ax.set_ylabel('mm')
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.25)
m.drawparallels(np.arange(-40.,-25.,5.),labels=[1,0,0,0])
m.drawmeridians(np.arange(-65.,-50.,5.),labels=[0,0,0,1],fontsize=10)
plt.ylabel('Latitud', fontsize=20)
plt.xlabel('Longuitud', fontsize=20)
plt.title('Precipitacion Estimada Acumulada \n Periodo 19-23/12/2015 \n por GPM IMERG (mm)', fontsize=25)


plt.savefig( 'PPestimadaPeriodo.png', dpi=300, format='png')



#########################################################################
#Acumulado Diario

# para abrir varios en una misma carpeta (por fecha):
#Fecha='20151223'
#DataPathd = '/data/miglesias/IMERGFR/'+Fecha+'/'
#Filelistd = sorted(os.listdir(DataPathd))

#ppestd=np.ma.zeros([len(Filelistd),152,152])	# [tiempos,lon,lat]
#for nfile in np.arange(0,len(Filelistd),1):
#	f = h5py.File(DataPathd + Filelistd[nfile], 'r')
#	a_group_key = list(f.keys())[0]
#	datad = f[a_group_key]
#	ppd = np.array(datad.get('precipitationCal'))
#	ppestd[nfile,:,:] = ppd[0,1149:1301,499:651] #aca estan en lon lat todavia
#	ppestd[ppestd < 0] =0

#ppacumd = np.sum(ppestd, axis=0)	# sumo todos los tiempos, acumulado total del periodo
#ppacumd = ppacumd*0.5   # multiplico por 0.5 para obtener mm/h
#ppacumd = np.transpose(ppacumd)	# necesito [lat, lon] y los datos vienen en [lon, lat]
# en el script de NCL de Cyn tambien intercambia el orden de lat, lon mediante pp = data(lat|:,lon|:)

#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap, cm

#plt.figure(figsize=(10,10))
#intervalos = [1,5,10,25,50,100,150,200,250]
#latcorners = ([-40,-25])
#loncorners = ([-65,-50])
#m = Basemap(projection='cyl',llcrnrlat=latcorners[0],urcrnrlat=latcorners[1],llcrnrlon=loncorners[0],urcrnrlon=loncorners[1],resolution='l')
#plt.contourf(lon,lat,ppacumd,intervalos,cmap=plt.cm.viridis_r,extend="both")
#cbar = plt.colorbar(ticks=intervalos,shrink=0.5)
#cbar.ax.set_ylabel('mm/hr')
#m.drawcoastlines(linewidth=0.5)
#m.drawcountries(linewidth=0.5)
#m.drawstates(linewidth=0.25)
#m.drawparallels(np.arange(-40.,-25.,5.),labels=[1,0,0,0])
#m.drawmeridians(np.arange(-65.,-50.,5.),labels=[0,0,0,1],fontsize=10)
#plt.ylabel('Latitud', fontsize=20)
#plt.xlabel('Longuitud', fontsize=20)
#Fecha=Fecha[0:4]+'_'+Fecha[4:6]+'_'+Fecha[6:8]
#plt.title('Precipitacion Estimada Acumulada Diario \n '+Fecha+'  \npor GPM IMERG (mm/hr)', fontsize=25)

#plt.savefig( 'PPestimadaDiaria'+ Fecha +'.png', dpi=300, format='png')




#########################################################################
#Acumulado cada 6 hrs
#A la fecha le introduzco un delta de tiempo de 6hs-> en este caso Se hizo juntando cada 6 horas previas los datos. Es decir tener carpetas de HORA 06 12 18 00  donde 

#00 06 tenia de 00.00-00.29.59 hasta 05.30.00-05.59.59
#en este caso 00 corresponde al acumulado de 18 a 00 del dia de la fecha en la que uno se para. 

# para abrir varios_:
#Fecha='20151224'
#Hora='00'
#DataPathh = '/data/miglesias/IMERGFR/HorarioAcu/'+Fecha+'/'+Hora+'/'
#Filelisth = sorted(os.listdir(DataPathh))
#ppesth=np.ma.zeros([len(Filelisth),152,152])	# [tiempos,lon,lat]
#for nfile in np.arange(0,len(Filelisth),1):
#	f = h5py.File(DataPathh + Filelisth[nfile], 'r')
#	a_group_key = list(f.keys())[0]
#	datah = f[a_group_key]
#	pph = np.array(datah.get('precipitationCal'))
#	ppesth[nfile,:,:] = pph[0,1149:1301,499:651] #aca estan en lon lat todavia
#	ppesth[ppesth < 0] =0

#ppacumh = np.sum(ppesth, axis=0)	# sumo todos los tiempos, acumulado total del periodo
#ppacumh = ppacumh*0.5   # multiplico por 0.5 para obtener mm/h
#ppacumh = np.transpose(ppacumh)	# necesito [lat, lon] y los datos vienen en [lon, lat]
# en el script de NCL de Cyn tambien intercambia el orden de lat, lon mediante pp = data(lat|:,lon|:)

#hay que guardar los datos cada 6 horas de pp acu para luego interpolarlos. hay que guardar en un npy los lat lon

#np.save('latIMERG.npy',lat)
#np.save('lonIMERG.npy',lon)

#np.ma.dump(ppacumh, '/data/miglesias/IMERGFR/DataImerg/'+Fecha+'_'+Hora)






















#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap, cm

#plt.figure(figsize=(10,10))
#intervalos = [1,5,10,25,50,100,150,200,250]
#latcorners = ([-40,-25])
#loncorners = ([-65,-50])
#m = Basemap(projection='cyl',llcrnrlat=latcorners[0],urcrnrlat=latcorners[1],llcrnrlon=loncorners[0],urcrnrlon=loncorners[1],resolution='l')
#plt.contourf(lon,lat,ppacumh,intervalos,cmap=plt.cm.YlGnBu,extend="both")
#cbar = plt.colorbar(ticks=intervalos,shrink=0.5)
#cbar.ax.set_ylabel('mm/hr')
#m.drawcoastlines(linewidth=0.5)
#m.drawcountries(linewidth=0.5)
#m.drawstates(linewidth=0.25)
#m.drawparallels(np.arange(-40.,-25.,5.),labels=[1,0,0,0])
#m.drawmeridians(np.arange(-65.,-50.,5.),labels=[0,0,0,1],fontsize=10)
#plt.ylabel('Latitud', fontsize=20)
#plt.xlabel('Longuitud', fontsize=20)
#Fecha=Fecha[0:4]+'_'+Fecha[4:6]+'_'+Fecha[6:8]
#if Hora=='00':
#	Hora='Desde 18-24hs '
#elif Hora=='06':
#	Hora='Desde 00-06hs '
#elif Hora=='12':
#	Hora='Desde 06-12hs '
#elif Hora=='18':
#	Hora='Desde 12-18hs '

#plt.title('Precipitacion Estimada Acumulada\n '+ Hora + Fecha +' \npor GPM IMERG (mm/hr)', fontsize=25)

#plt.savefig( 'PPestimadaHoraria'+ Hora + Fecha +'.png', dpi=300, format='png')



