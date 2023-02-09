#
#	m-a_sep_TNO.py
#
#	Makes a 2d histogram of the mass semi-major axis separation distribution given a list of J2 vals
#	with data input specifically from Borasisi. This could in principle be replaced with other J2 values
#	Additional changes to indicate hierarchical triple detectability and for input of any TNO system by William Giforos, 09/2022
#
#	Benjamin Proudfoot 
#	06/13/22
#

from unpack import *
import numpy as np
import matplotlib.pyplot as plt
import corner
import sys
import matplotlib.colors as colors
import os
import pandas as pd
import commentjson as json
import astropy.io.fits

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

# Ensure that its being run from the right place
if 'results' in os.getcwd():
    runprops = ReadJson('runprops.txt').outProps()
    objectdata = ReadJson('objectdata.txt').outProps()
    resultspath = os.getcwd()
    print(resultspath)
else:
    print("m-a_sep_borasisi.py needs to be run from a results directory.")
    sys.exit()

name_dict = runprops.get("names_dict")
objectnames = []
for i in name_dict.values():
    objectnames.append(i)
print(objectnames)

# Import multimoon sigsdf file
sigsdf_mm = pd.read_csv(resultspath +'/sigsdf_mm.csv',sep=',',index_col=0)
#print(sigsdf_mm)

# Define the total mass, mass ratio, and fractional brightness (for purposes of the Stability Limit)
Mtot = (sigsdf_mm.at['mass_1','median'] + sigsdf_mm.at['mass_2','median']) * 10**(18) # or 3.43 E18 kg (borasisi & pabu)
# Eris data: Total system mass 1.66 *10^22
# Eris proposed mass ratio .05, sma 12,400 km
#Mtot = 1.5*(10**22)
print("M_tot:",Mtot," kg")
x =  .05# mass ratio, m in / m primary
fr_br = 10**(objectdata.get('dmag')*0.4) # Convert dmag to fractional brightness

# solve for the masses
Mout = Mtot / (((1+x) / np.sqrt((fr_br*(1 + x**(2/3)))**3)) + 1)
print("M_out:",Mout," kg")
Mp = (Mtot - Mout) / (1 + x)
print("M_prim:",Mp," kg")
Min = x*Mp
print("M_in:",Min," kg")

# Calculate the inner semimajor axis using the j2r2 median value
a_out = objectdata.get('sma_outer') #km
a_in = np.sqrt(2*sigsdf_mm.at['j2r2_1','median']*np.reciprocal(x))
# Eris semimajor axis
#a_in = 12400
print("a_in:",a_in," km")

# Calculate the period of the inner orbit in days
G = 6.6742867 * 10**(-11) # N*(m^2)(kg^-2)
P = np.sqrt((4*((np.pi)**2)*((a_in*1000)**3)) / (G*(Mp+Min)))
P = P*(1/60)*(1/60)*(1/24) # convert to days
print("P_in:",P," days")

# Define the distance that 1 WFC3 pixel covers
ang_size = 40 / 1000 / 3600 # from mas to deg
pix_dist_1 = np.tan(np.deg2rad(ang_size)) * objectdata.get('avg_distance') * (1.496 * (10**8)) # converts AU to km as well
print("WFC3 pixel:", pix_dist_1," km/pixel")

# Calculate for 3 pixels
ang_sz_3 = ang_size * 3
pix_dist_3 = np.tan(np.deg2rad(ang_sz_3)) * objectdata.get('avg_distance') * (1.496 * (10**8)) # converts AU to km as well

# Define the separation between the inner moon and the primary object in units of pixels
sep_pix = a_in / pix_dist_1 
print("sep_pix_inner:", sep_pix," pixels")

# Inputting j_2*R^2 values
signeg1 = sigsdf_mm.at['j2r2_1','median'] + sigsdf_mm.at['j2r2_1','-1sigma']
signeg2 = sigsdf_mm.at['j2r2_1','median'] + sigsdf_mm.at['j2r2_1','-2sigma']
signeg3 = sigsdf_mm.at['j2r2_1','median'] + sigsdf_mm.at['j2r2_1','-3sigma']
sig1 = sigsdf_mm.at['j2r2_1','median'] + sigsdf_mm.at['j2r2_1','1sigma']
sig2 = sigsdf_mm.at['j2r2_1','median'] + sigsdf_mm.at['j2r2_1','2sigma']
sig3 = sigsdf_mm.at['j2r2_1','median'] + sigsdf_mm.at['j2r2_1','3sigma']
med = sigsdf_mm.at['j2r2_1','median']
j2r2sigs = [signeg3, signeg2, signeg1, med, sig1, sig2, sig3]
#j2r2sigs = [1374.529928,2508.523643,4689.973173,7705.891633,12105.60188,17474.50743,21241.29693]
#j2r2sigs = [1374.529928,2508.523643,4689.973173,12105.60188,17474.50743,21241.29693]	# List without median included
j2r2sigs = np.array(j2r2sigs)

# Making latex compatible labels for nice looking figures
labels = [r"$-$3$\sigma$",r"$-$2$\sigma$",r"$-$1$\sigma$",r"Median",r"1$\sigma$",r"2$\sigma$",r"3$\sigma$"]

# Define the graph ceiling value
if a_out < 10000:
    y_max = 10000
    mrl_vertical = np.log10((y_max - pix_dist_3)*.001)
    label_shift = 200
else:
    y_max = 100000
    mrl_vertical = .59
    label_shift = 800

# Setting up arrays for mass ratio and sma and filling a grid of J2 values
m2m1 = np.linspace(0,1,num = 1000) #x-axis
sma = np.linspace(1,y_max,num = 1000) #y-axis
grid = np.empty((1000,1000))
for i,mr in enumerate(m2m1):
	grid[i,:] = 0.5*(mr/(1+mr)**2)*(sma**2)
j2r2s = grid.transpose()

# Load in fonts (will not work without my specific font set up)
from matplotlib import font_manager
#font_path = "/Users/benjamin.proudfoot/Documents/Research/.mm_venv/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/HelveticaNeueLight.ttf"
#font_path = "/home/byu.local/wgiforos/.vnpsf/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/HelveticaNeueLight.ttf"
#font_manager.fontManager.addfont(font_path)
#prop = font_manager.FontProperties(fname=font_path)
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = prop.get_name()

# Filling in colored contours
fig,ax = plt.subplots()
CS = ax.contourf(m2m1,sma,j2r2s,j2r2sigs, colors = ["#481567FF","#20A387FF","#FDE725FF","#FDE725FF","#20A387FF","#481567FF"], alpha = 0.8)

# Filling in median J2 line
ax.contour(m2m1,sma,j2r2s,[sigsdf_mm.at['j2r2_1','median']], colors = "black", zorder = 1)

# Setting up lines off the plot so I can have a legend the way I want it
dummyx = np.linspace(0,1, num = 10)
dummyy = -1000*dummyx
plt.plot(dummyx, dummyy, color = "black", label = "Median $J_2$")
plt.plot(dummyx, dummyy, color = "#FDE725FF", label = r"1$\sigma$ limit", linewidth = 5)
plt.plot(dummyx, dummyy, color = "#20A387FF", label = r"2$\sigma$ limit", linewidth = 5)
plt.plot(dummyx, dummyy, color = "#481567FF", label = r"3$\sigma$ limit", linewidth = 5)
plt.legend(loc = "lower right")

# Fiddling with plotting settings
plt.xlabel(r"$\frac{M_{in}}{M_{prim}}$", fontsize = 20)
plt.ylabel(r"$a$ (km)", fontsize = 20)
plt.title(objectnames[0] + " 3rd Body Prediction")
# Eris plot x axis
#plt.xlim(0.0000001,1)
plt.xlim(0.0001,1)
plt.ylim(40,y_max)
plt.yscale("log")
plt.xscale("log")
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.gca().invert_xaxis()

# Not sure what this does.... I forget
fmt = {}
for l,s in zip(CS.levels, labels):
	fmt[l] = s

# Plot points for comparison bodies
m1 = 7.657
m1e = (0.323 + 0.346)/2
m2 = 5.725
m2e = (0.244 + 0.234)/2
a = 838
aep = 13
aem = 21

m2m1lempo = m2/m1
m2m1lempo_err = (m2/m1)*np.sqrt( (m2e/m2)**2 + (m1e/m1)**2 )

plt.scatter(m2m1lempo, 838, s = 15, color = "black")
plt.text(0.7, 650, "Lempo-Hiisi", rotation = 0)

plt.scatter((138/181)**3, 349, s = 15, color = "black")
plt.text(0.93, 380, "WC510", rotation = 0)

plt.scatter(1.0, 198.6, s = 15, color = "black")
plt.text(0.9, 170, "JY31", rotation = 0)

# Calculate periapsis
a = objectdata.get('sma_outer')
e = objectdata.get('ecc_outer')
periapsis = a*(1-e) 

"""
# Brightness Ratio limit (from which we calculate a mass ratio limit)
tot_counts = objectdata.get('total_counts')
sky_noise = objectdata.get('sky_noise')
f = "../../" + runprops.get('image_path') + runprops.get('input_image')
fitsfile = astropy.io.fits.open(f)
SNR = fitsfile[1].header['SNRMEAN']
br = (sky_noise * SNR) / tot_counts
mr_lim = br**(3/2)
print("mr_lim:",mr_lim)

# Plot Brightness Ratio limit (Mass Ratio limit)
plt.axvline(x = mr_lim, ymin = np.log10((y_max - pix_dist_3)*.001), ymax = 1, linestyle = "--", color = "black", alpha = 0.6, zorder = 1)
#plt.text(0.001, 650, "brightness")
#plt.text(0.001, 500, "detection")
#plt.text(0.001, 380, "limit")
"""
# Brightness Ratio limit (from which we calculate a mass ratio limit) using 10 instead of SNR
tot_counts = objectdata.get('total_counts')
sky_noise = objectdata.get('sky_noise')
br = (sky_noise * 10) / tot_counts
mr_lim_10 = br**(3/2)
print("mr_lim_10:",mr_lim_10)

# Plot Brightness Ratio limit (Mass Ratio limit)
plt.axvline(x = mr_lim_10, ymin = mrl_vertical, ymax = 1, linestyle = "--", color = "black", alpha = 0.6, zorder = 1)
plt.text(mr_lim_10 - .0001,  pix_dist_3 + 5000, "approximate")
plt.text(mr_lim_10 - .0001,  pix_dist_3 + 2300, "detectability")
plt.text(mr_lim_10 - .0001,  pix_dist_3 + 500, "limit")

# Eris text line
#plt.axvline(x = mr_lim_10, ymin = .67, ymax = 1, linestyle = "--", color = "black", alpha = 0.6, zorder = 1)
#plt.text(mr_lim_10 - .0000001,  pix_dist_3 + 7000, "approximate")
#plt.text(mr_lim_10 - .0000001,  pix_dist_3 + 3000, "detectability")
#plt.text(mr_lim_10 - .0000001,  pix_dist_3 + 400, "limit")

# plot beta object periapsis
plt.axhline(y = periapsis, xmin = 0.0001, xmax = 1.0, linestyle = "--", color = "black", alpha = 0.4, zorder = 1)
plt.text(0.9, periapsis + label_shift, objectnames[1] + " periapsis")

# plot 1 WFC3 pixel distance
plt.axhline(y = pix_dist_1, xmin = 0.0001, xmax = 1.0, linestyle = "--", color = "black", alpha = 0.4, zorder = 1)
plt.text(0.9, pix_dist_1 + 100, r"$\sim$1 WFC3 pixel")

# plot 3 WFC3 pixel distance
plt.axhline(y = pix_dist_3, xmin = 0.0001, xmax = 1.0, linestyle = "--", color = "black", alpha = 0.4, zorder = 1)
plt.text(0.9, pix_dist_3 + 200, r"$\sim$3 WFC3 pixel")

# plot primary radius
plt.axhline(y = objectdata.get('primary_radius'), xmin = 0.0001, xmax = 1.0, linestyle = "--", color = "black", alpha = 0.4, zorder = 1)
plt.text(0.026, objectdata.get('primary_radius') + 5, objectnames[0] + " radius")

# Shade apropriate regions
plt.fill_between(x = np.linspace(0,1,num = 1000), y1 =  y_max, y2 = periapsis, alpha = 0.5, color = 'dimgray')

plt.fill_betweenx(y = np.linspace(1,y_max,num = 1000), x1 =  mr_lim_10, x2 = 1, where = np.linspace(1,y_max,num = 1000) >= pix_dist_3, alpha = 0.5, color = 'dimgray')

# Saving figure
plt.tight_layout()
plt.savefig("mr_x_sma_" + objectnames[0] + ".pdf")
