import numpy as np
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
                                     ra_units='degrees', dec_units='degrees', npatch=50)
        self.rand = treecorr.Catalog(ra=randoms['ra'], dec=randoms['dec'], 
                                     ra_units='degrees', dec_units='degrees', npatch=50)

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
            
            
    