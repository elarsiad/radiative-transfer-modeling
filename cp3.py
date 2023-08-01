#Computational Problem 2
import math
import numpy as np
import matplotlib.pyplot as plt

### Computing albedo spectrum
## Read the albedo files
f = open(r'data\albedo.txt','r')
lines = f.readlines()
#a. Wavelength (nm)
#b. Surface reflectance values for sand surface
#c. Surface reflectance values for soil surface
#d. Surface reflectance values for snow surface
#e. Surface reflectance values for vegetation surface
#f. Surface reflectance values for water surface
#
#     a       b       c       d       e       f
# 280.0   0.091   0.030   0.949   0.059   0.025
del lines[0:8]

wls = []
alb_sand = []
alb_soil = []
alb_snow = []
alb_vegt = []
alb_watr = []

for line in lines:
    wls.append(float(line[0:6]))
    alb_sand.append(float(line[9:14]))
    alb_soil.append(float(line[17:22]))
    alb_snow.append(float(line[25:30]))
    alb_vegt.append(float(line[33:38]))
    alb_watr.append(float(line[41:46]))
f.close()

alb_all = np.transpose(np.array([alb_sand,alb_soil,alb_snow,alb_vegt,alb_watr]))
alb_mean = ( 5.8*alb_all[:,0] + 
             7.5*alb_all[:,1] + 
             5.8*alb_all[:,2] +
             9.1*alb_all[:,3] + 
            70.8*alb_all[:,4] )/99

plt.figure(1)
plt.plot(wls,alb_watr)
plt.plot(wls,alb_sand)
plt.plot(wls,alb_vegt)
plt.plot(wls,alb_soil)
plt.plot(wls,alb_snow)
plt.plot(wls,alb_mean,'--')
plt.legend([r'water (70.8%)',
            r'sand (5.8%)',
            r'vegetation (9.1%)',            
            r'soil (7.5%)',
            r'snow (5.8%)',
            r"earth's average"])
plt.ylabel('Albedo')
plt.xlabel('Wavelength [nm]')
plt.grid()

## Read the atmosphere layers effective ztp from problem 1 of CP1
f = open('data\ztp_eff.txt','r')
lines = f.readlines()
f.close()

z = []
t = []
p = []

nlayer = 0
for line in lines:
    z.append(float(line[0:4])) #km
    t.append(float(line[5:11])) #K
    p.append(float(line[12:23])) #hPa
    nlayer += 1

## Read mixing ratio from conc.txt
f = open('data\conc.txt','r')
lines = f.readlines()
f.close()
del lines[0] #remove header line

conc_dat = []
for line in lines:
    linetemp = [line[0:5],line[5:17],line[17:29],line[29:41],line[41:53],line[53:65],line[65:77]]
    conc_dat.append([float(i) for i in linetemp])

## Read the absorption cross section of H2O, O2, and O3
siga_h2o = []
siga_o2 = []
siga_o3 = []

f = open('data\h2o.txt','r')
lines = f.readlines()
for line in lines:
    siga_h2o.append(float(line[8:18]))
f.close()

f = open('data\o2.txt','r')
lines = f.readlines()
for line in lines:
    siga_o2.append(float(line[8:18]))
f.close()

f = open('data\o3.txt','r')
lines = f.readlines()
for line in lines:
    siga_o3.append(float(line[8:18]))
f.close()

SZA = 30
#wls = range(280,1001) # (nm) wavelength

mr = [[0 for i in range(nlayer)] for j in wls] # shape: 721 x 10
incang = [[0 for i in range(nlayer)] for j in wls] # shape: 721 x 10

## Computing incident angle and path length per layer for each wavelength
PL = [[0 for i in range(nlayer)] for j in wls] # shape: 721 x 10
PLsum = [0 for i in wls] # shape: 721

for j,wl in enumerate(wls):
    wl = wl * 1e-3 # convert wavelength (nm to micron)
    for i in reversed(range(nlayer)):
        mr[j][i] = (6432.8 + 2949810/(146-wl**(-2)) + 25540/(41-wl**(-2))) / (1013.25/288.15) * (p[i]/t[i]) * 1e-8 + 1
        if i == 9:
            incang[j][i] = SZA
        else:
            incang[j][i] = np.arcsin(mr[j][i+1]*math.sin(incang[j][i+1]/180*math.pi)/mr[j][i]) *180/math.pi
        PL[j][i] = 5/math.cos(incang[j][i]/180*math.pi)
        PLsum[j] = PLsum[j] + PL[j][i]
    PLsum[j] = 50/math.cos(SZA/180*math.pi) - PLsum[j]

## Computing radiance

R  = 287.058;   #(J/(kg K)) Specific gas constant for air
Av = 6.02295 * 1e23 #(/mol) Avogadro number
M = 28.97e-3 #(kg/mol) molecular weight of air
delt = 0.035
fdelt = (6+3*delt)/(6-7*delt)
#A = 1.0 #Albedo

def fracrad(flag_scat,flag_abs,flag_cloud,A):

    UWbot = [0 for i in wls] # shape: 721
    UWsum = [0 for i in wls] # shape: 721
    UW = [[0 for i in range(nlayer)] for j in wls] # shape: 721 x 10
    DW = [[0 for i in range(nlayer)] for j in wls] # shape: 721 x 10

    def light_scattered(I_in):
        sigs = 8*math.pi**3*(mr[j][i]**2-1)**2/(3*wl**4*Ns**2)*fdelt
        I_scat = I_in*(1-math.exp(-sigs*Ns*PL[j][i]*1e5))
        I_scatfor = I_scat*3/4*(1+math.cos(0/180*math.pi)**2)/(4*math.pi)
        return (I_scat,I_scatfor)

    def light_absorbed(I_in):
        Cx = conc_dat[int(z[i])*10][1] #h2o conc (v/v)
        Na_h2o = Cx*Ns #(/cm3) number density (n/v)
        Cx = conc_dat[int(z[i])*10][2] #o3 conc (v/v)
        Na_o3 = Cx*Ns
        Na_o2 = 0.21*Ns
        I_abs = 0
        for siga,Na in [(siga_h2o[j], Na_h2o), (siga_o2[j],Na_o2), (siga_o3[j],Na_o3)]:
        #for siga,Na in [(siga_o2[j],Na_o2)]:
            I_abs = I_abs + I_in*(1-math.exp(-siga*Na*PL[j][i]*1e5))
        return I_abs
        
    alb_cloud = 0.8

    for j,wl in enumerate(wls):
        wl = wl * 1e-7 # convert wavelength (nm to cm)
        for i in reversed(range(nlayer)):
            rho = p[i]*100/t[i]/R/1000 #(g/cm3)
            Ns = rho*Av/M*1e-3 #(/cm3)
            if i == 9:
                Iin = 1*math.cos(SZA/180*math.pi)
            elif i == 0 and flag_cloud is True:
                Iin = DW[j][i+1]*math.cos(incang[j][i+1]/180*math.pi)/math.cos(incang[j][i]/180*math.pi)
                Iin = Iin*0.7
                UW[j][i] = UW[j][i] + Iin/0.7*0.3*alb_cloud/math.pi
            else:
                Iin = DW[j][i+1]*math.cos(incang[j][i+1]/180*math.pi)/math.cos(incang[j][i]/180*math.pi)
            #scattered light
            Iscat,Iscatfor = light_scattered(Iin)
            #source term SZA of the scattered light as upwelling
            UW[j][i] = UW[j][i] + Iscat*3/4*(1+math.cos((180-incang[j][i])/180*math.pi)**2)/(4*math.pi)
            #absorbed light
            Iabs = light_absorbed(Iin)
            #Incoming - scattered + forward scattered as downwelling - absorbed
            if flag_scat is False:
                Iscat = 0
                Iscatfor = 0
                UW[j][i] = 0
            if flag_abs is False:
                Iabs = 0
            DW[j][i] = Iin - Iscat + Iscatfor - Iabs
            if DW[j][i] < 0:
                DW[j][i] = 0
        UWbot[j] = DW[j][i]*A[j]/math.pi

    for j,wl in enumerate(wls):
        wl = wl * 1e-7 # convert wavelength (nm to cm)
        UWsum[j] = UWbot[j]
        for i in range(nlayer): 
            rho = p[i]*100/t[i]/R/1000 #(g/cm3)
            Ns = rho*Av/M*1e-3 #(cm3)
            #scattered light
            Iscat,Iscatfor = light_scattered(UWsum[j])
            #absorbed light
            Iabs = light_absorbed(UWsum[j])
            #Incoming UW - scattered + forward scattered - absorbed + UW each layer
            if flag_scat is False:
                Iscat = 0
                Iscatfor = 0
            if flag_abs is False:
                Iabs = 0
            UWsum[j] = UWsum[j] - Iscat + Iscatfor - Iabs
            #if i == 0 and flag_cloud is True:
            #    UWsum[j] = UWsum[j]*0.7
            UWsum[j] = UWsum[j] + UW[j][i]
            if UWsum[j] < 0:
                UWsum[j] = 0
                
    return [iUWsum*np.pi/np.cos(SZA/180*np.pi) for iUWsum in UWsum]

UWsum_both = fracrad(True,True,False,alb_mean)
UWsum_scat = fracrad(True,False,False,alb_mean)
UWsum_abs = fracrad(False,True,False,alb_mean)

## Read solar spectrum from solar.txt
f = open('data\solar.txt','r')
lines = f.readlines()

solar_dat = []
for line in lines:
    linetemp = float(line[7:18])
    solar_dat.append(linetemp)
f.close()

#plot solar spectrum
plt.figure(2)
plt.plot(wls,solar_dat)
plt.title('Solar spectrum')
plt.xlabel('wavelength [nm]')
plt.ylabel('[W/m2/nm]')
plt.xlim([280,1000])
plt.grid()

#plot albedo spectrum
plt.figure(3)
plt.title('Albedo spectrum')
plt.plot(wls,UWsum_both)
plt.plot(wls,UWsum_scat)
plt.plot(wls,UWsum_abs)
plt.legend(['scattering and absorption','pure scattering','pure absorption'])
plt.xlabel('wavelength[nm]')
plt.ylabel('Albedo')
plt.xlim([280,1000])
plt.grid()

#plot solar radiance spectrum
albspec_both = [0 for i in wls]
albspec_scat = [0 for i in wls]
albspec_abs = [0 for i in wls]
for j,wl in enumerate(wls):
    albspec_both[j] = UWsum_both[j]*solar_dat[j]
    albspec_scat[j] = UWsum_scat[j]*solar_dat[j]
    albspec_abs[j] = UWsum_abs[j]*solar_dat[j]


#plot spectral resolution effect
plt.figure(4)
plt.title('Spectral resolution effect on the received solar radiance spectrum')
for i,specr in enumerate([1,10,100]):
    plt.subplot(311+i)
    plt.plot(wls[::specr],albspec_both[::specr])
    plt.plot(wls[::specr],albspec_scat[::specr])
    plt.plot(wls[::specr],albspec_abs[::specr])
    plt.text(800,0.25,'Resolution '+str(specr)+' nm')
    plt.ylabel('[W/m2/nm]')
    plt.xlim([280,1000])
    plt.grid()
plt.legend(['scattering and absorption','pure scattering','pure absorption'])


#plot cloud effect
UWsum_both_cloud = fracrad(True,True,True,alb_mean)
albspec_both_cloud = [0 for i in wls]
for j,wl in enumerate(wls):
    albspec_both_cloud[j] = UWsum_both_cloud[j]*solar_dat[j]

plt.figure(5)
plt.title('Albedo spectrum')
plt.plot(wls,UWsum_both)
plt.plot(wls,UWsum_both_cloud)
plt.legend(['not cloudy','30% cloudy'])
plt.xlabel('Wavelength [nm]')
plt.ylabel('Albedo')
plt.xlim([280,1000])
plt.grid()


#plot surface effect
UWsum_both_snow = fracrad(True,True,False,alb_snow)
UWsum_both_vegt = fracrad(True,True,False,alb_vegt)
UWsum_both_watr = fracrad(True,True,False,alb_watr)
albspec_both_snow = [0 for i in wls]
albspec_both_vegt = [0 for i in wls]
albspec_both_watr = [0 for i in wls]
for j,wl in enumerate(wls):
    albspec_both_snow[j] = UWsum_both_snow[j]*solar_dat[j]
    albspec_both_vegt[j] = UWsum_both_vegt[j]*solar_dat[j]
    albspec_both_watr[j] = UWsum_both_watr[j]*solar_dat[j]

plt.figure(6)
plt.title('Albedo spectrum')
plt.plot(wls,UWsum_both,'--')
plt.plot(wls,UWsum_both_snow)
plt.plot(wls,UWsum_both_vegt)
plt.plot(wls,UWsum_both_watr)
plt.legend(['average earth','snow','vegetation','ocean'])
plt.xlabel('Wavelength [nm]')
plt.ylabel('Albedo')
plt.xlim([280,1000])
plt.grid()
plt.show()
