import astropy.io.fits as fits
import numpy as np
from scipy.integrate import quad
import scipy.integrate as integrate
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table,join
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.cosmology import FlatLambdaCDM
import treecorr
#import arviz as az
#import pymc as pm

#%%
home_dir = os.path.expanduser('~')
# Construct the path to the "Thesis" directory on the desktop
thesis_path = os.path.join(home_dir, 'Desktop', 'Thesis')


fits_file_path = os.path.join(thesis_path, "table")  


t= Table.read(fits_file_path)





fits_random = os.path.join(thesis_path, "SN-C3_randoms_ugriz_trim_video.fits") 

hdulist = fits.open(fits_random)
hdulist.info()

t2= Table.read(fits_random)




config = {
    'min_sep':  0.003,
    'max_sep': 1.78,
    'bin_size': 0.1,
    'sep_units': 'degrees',
    'var_method': 'bootstrap'  # or 'jackknife'
}


#%%

class Subsample:
    def __init__(self, z_min, z_max, SM_min, SM_max):
        self.z_min = z_min
        self.z_max = z_max
        self.SM_min = SM_min
        self.SM_max = SM_max

    def apply(self, data):
        """
        Selects data points within the subsample region.

        Args:
            data: A dictionary containing data columns like 'z' and 'SM'.

        Returns:
            A boolean mask for the subsample selection.
        """
        return (data['z'] > self.z_min) & (data['z'] <= self.z_max) & \
               (data['SM'] > self.SM_min) & (data['SM'] <= self.SM_max)
    
    
    
    

# Define z and SM ranges
z_values = np.array([ 0.8, 1.0, 1.2, 1.4])
SM_range = np.linspace(9, 11, num=3)

# Create subsamples
subsamples = []
z_mean_range = []
SM_mean_range = []

for i in range(len(z_values) - 1):
    z_min = z_values[i]
    z_max = z_values[i + 1]
    z_mean = (z_min + z_max) / 2
    
    z_mean_range.append(z_mean)

    for j in range(len(SM_range) - 1):
        SM_min = SM_range[j]
        SM_max = SM_range[j + 1]
        SM_mean = (SM_min + SM_max) / 2
        
        SM_mean_range.append(SM_mean)

        
        subsample = Subsample(z_min, z_max, SM_min, SM_max)
        subsamples.append(subsample)

        
        
plt.scatter(t['z'], t['SM'], label='All galaxies')

# Scatter plot with subsample colors
colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown', 'pink', 'gray', 'yellow', 'cyan']

for i, subsample in enumerate(subsamples):
    subset = subsample.apply(t)
    z_subsample = t['z'][subset]
    SM_subsample = t['SM'][subset]

    plt.scatter(z_subsample, SM_subsample, label=f'Subsample {i+1}', c=colors[i])

for i in range(len(z_values) - 1):
    z_min = z_values[i]
    z_max = z_values[i + 1]
    plt.axvline(z_min, linestyle='--', color='gray')
    plt.axvline(z_max, linestyle='--', color='gray')

plt.xlabel('Redshift (z)')
plt.ylabel('Stellar Mass (SM)')
plt.title('Galaxy Mass Distribution')
plt.legend()
plt.show()

# Create SkyCoord catalogs for subsamples
catalogs = []

for subsample in subsamples:
    subset = subsample.apply(t)
    catalog = t[subset]

    #catalog = SkyCoord(ra=ra_subset * u.deg, dec=dec_subset * u.deg)
    catalogs.append(catalog)

# Print number of galaxies in each subsample
for i, catalog in enumerate(catalogs):
    N = len(catalog)
    print(f"Subsample {i+1}: N = {N}")



class CorrelationFunction:
    def __init__(self, cat, randoms, config, rr=None):
        self.cat = cat
        self.randoms = randoms
        self.config = config

        # Create TreeCorr catalogs
        self.data = treecorr.Catalog(ra=cat['ra'], dec=cat['dec'], 
                                     ra_units='degrees', dec_units='degrees', npatch=50)
        self.rand = treecorr.Catalog(ra=randoms['ra'], dec=randoms['dec'], 
                                     ra_units='degrees', dec_units='degrees', npatch=50)

        # Create TreeCorr correlation functions
        self.dd = treecorr.NNCorrelation(**config)
        self.rr = treecorr.NNCorrelation(**config)
        self.dr = treecorr.NNCorrelation(**config)
        
    def process_rand(self):
        self.rr.process(self.rand) ###????maybe here problem?

    def process(self):
        print("a")
        self.dd.process(self.data)
        print("b")
        self.dr.process(self.data, self.rand)

    def calculate_w_theta(self):
        theta = np.exp(self.dd.meanlogr)
        w, varw = self.dd.calculateXi(rr=self.rr, dr=self.dr) 
        
        return w, varw, theta , self.rr
    

    def power_law(self,theta, A, delta=-0.8):
        self.A=A
        self.theta=theta
        
        return A * theta**delta


#%%


subsample = (t['z'] > 0.5) & (t['z'] <= 0.6) & (t['SM'] > 10) & (t['SM'] <= 10.5)
subset = t[subsample]
randoms = t2[::1000]

corr_func = CorrelationFunction(subset, randoms, config)



#%%


corr_func.process_rand() #crashes here
corr_func.process() # and crashes here too
    
#np.save('corr_func.npy', corr_func.dd, corr_func.dr, corr_func.rr)


np.save('corr_func_dd.npy', corr_func.dd)
np.save('corr_func_dr.npy', corr_func.dr)
np.save('corr_func_rr.npy', corr_func.rr)
