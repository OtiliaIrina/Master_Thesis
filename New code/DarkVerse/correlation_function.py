import numpy as np
import treecorr

class CorrelationFunction:
    def __init__(self, cat, randoms, config, rr=None):
        self.cat = cat
        self.randoms = randoms
        self.config = config
        self.rr_precomputed = rr  # Store precomputed RR pairs 

        # Create TreeCorr catalogs
        self.data = treecorr.Catalog(ra=cat['ra'], dec=cat['dec'], 
                                     ra_units='degrees', dec_units='degrees', npatch= 1)
        self.rand = treecorr.Catalog(ra=randoms['ra'], dec=randoms['dec'], 
                                     ra_units='degrees', dec_units='degrees', npatch= 1)

        # Create TreeCorr correlation functions
        self.dd = treecorr.NNCorrelation(**config)
        self.dr = treecorr.NNCorrelation(**config)
        self.rr = treecorr.NNCorrelation(**config) if rr is None else rr  # Use precomputed RR if provided

    def process(self):
        self.dd.process(self.data)
        self.dr.process(self.data, self.rand)
        
        if self.rr_precomputed is None:  # Compute RR only if not provided
            self.rr.process(self.rand)
            self.rr_precomputed = self.rr  # Store computed RR to avoid reprocessing

    def calculate_w_theta(self):
        """Calculate the angular correlation function."""
        theta = np.exp(self.dd.meanlogr)
        w, varw = self.dd.calculateXi(rr=self.rr_precomputed, dr=self.dr) 
        
        return w, varw, theta, self.rr_precomputed, self.dr, self.dd

    def bootstrap_w_theta(self, num_bootstrap=100):
        """Perform bootstrap resampling and compute variance/covariance of w_theta."""
        bootstrap_w_theta = []

        for _ in range(num_bootstrap):
            # Sample galaxy indices with replacement
            sample_indices = np.random.choice(len(self.cat), size=len(self.cat), replace=True)
            bs_galaxies = self.cat[sample_indices]

            # Create new TreeCorr catalog for resampled galaxies
            bs_data = treecorr.Catalog(ra=bs_galaxies['ra'], dec=bs_galaxies['dec'], 
                                       ra_units='degrees', dec_units='degrees', npatch=50)

            # Compute DD and DR for resampled data
            bs_dd = treecorr.NNCorrelation(**self.config)
            bs_dr = treecorr.NNCorrelation(**self.config)

            bs_dd.process(bs_data)
            bs_dr.process(bs_data, self.rand)

            # Calculate w_theta using precomputed RR
            w_theta_bs, varw_bs = bs_dd.calculateXi(rr=self.rr_precomputed, dr=bs_dr) 

            bootstrap_w_theta.append(w_theta_bs)

        # Convert list to NumPy array
        bootstrap_w_theta = np.array(bootstrap_w_theta)

        # Compute variance and covariance
        variance_bootstrap = np.var(bootstrap_w_theta, axis=0)
        covariance_bootstrap = np.cov(bootstrap_w_theta, rowvar=False)

        return variance_bootstrap, covariance_bootstrap

    def write_results(self, file_name):
        """Save the correlation function to a file."""
        self.dd.write(file_name, rr=self.rr_precomputed, dr=self.dr)

    def read_results(self, file_name, file_type=None):
        """Load correlation function results from a file."""
        self.dd.read(file_name, file_type=None)
        self.dr.read(file_name, file_type=None)
        self.rr.read(file_name, file_type=None)

    def print_num_pairs(self):
        """Print the number of galaxy pairs used in the correlation function."""
        print("Number of pairs (per bin):")
        print(f"DD pairs: {self.dd.npairs}")
        print(f"DR pairs: {self.dr.npairs}")
        print(f"RR pairs: {self.rr_precomputed.npairs}")

        
        
"""import numpy as np
from scipy.integrate import quad
import scipy.integrate as integrate
from astropy import units as u
from astropy.table import Table,join
import treecorr
import halomod as hm
import hmf





class CorrelationFunction:
    def __init__(self, cat, randoms, config, rr=None):
        self.cat = cat
        self.randoms = randoms
        self.config = config

        # Create TreeCorr catalogs
        self.data = treecorr.Catalog(ra=cat['ra'], dec=cat['dec'], 
                                     ra_units='degrees', dec_units='degrees', npatch = len(cat) // 2)
        self.rand = treecorr.Catalog(ra=randoms['ra'], dec=randoms['dec'], 
                                     ra_units='degrees', dec_units='degrees', npatch = len(cat) // 2)

        # Create TreeCorr correlation functions
        self.dd = treecorr.NNCorrelation(**config)
        self.rr = treecorr.NNCorrelation(**config)
        self.dr = treecorr.NNCorrelation(**config)
        

    def process(self):
        self.dd.process(self.data)
        self.dr.process(self.data, self.rand)
        self.rr.process(self.rand)


    def calculate_w_theta(self):
        theta = np.exp(self.dd.meanlogr)
        w, varw = self.dd.calculateXi(rr=self.rr, dr=self.dr) 
        
        return w, varw, theta ,self.rr, self.dr, self.dd
    
    def write_results(self, file_name):
            self.dd.write(file_name, rr=self.rr, dr=self.dr)
            
            
    def read_results(self, file_name, file_type=None):

        # Read correlation function data into the `dd`, `dr`, and `rr` objects
        self.dd.read(file_name, file_type=None)
        self.dr.read(file_name, file_type=None)
        self.rr.read(file_name, file_type=None)
        
        
    def print_num_pairs(self):
    
        print("Number of pairs (per bin):")
        print(f"DD pairs: {self.dd.npairs}")
        print(f"DR pairs: {self.dr.npairs}")
        print(f"RR pairs: {self.rr.npairs}")
            
            
    """