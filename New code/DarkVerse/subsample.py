from .correlation_function import CorrelationFunction
import numpy as np
from scipy.optimize import curve_fit
import halomod as hm

class Subsample:
    def __init__(self, catalog, randoms, z_min, z_max, SM_min, SM_max, config, w_theta=None, theta=None):
        self.catalog = catalog
        self.randoms = randoms
        self.z_min = z_min
        self.z_max = z_max
        self.SM_min = SM_min
        self.SM_max = SM_max
        self.config = config  # Configuration for correlation function
        self.info = {}  # Empty dictionary to store information

        # Compute the redshift distribution
        self.nz = hm.integrate_corr.flat_z_dist(self.z_min, self.z_max)

        self.mask = self.apply(catalog)  # Selection of the subsample
        self.filtered_catalog = catalog[self.mask]  # Subsample of galaxies
        self.N = len(self.filtered_catalog)  # Number of galaxies in the subsample
        
        # Initialize attributes for correlation function
        self.corr_func = None
        self.power_law_params = None  # (A, gamma) from power-law fit
        self.gg = None  # Placeholder for computed correlation functions using halomod(gg.angular_corr_matter and gg.angular_corr_gal)

        # Measure correlation function
        if w_theta is None or theta is None:
            self.measure_w_theta()
        else:
            self.w_theta = w_theta
            self.theta = theta

        # Compute the angular correlation functions
        self.compute_gg()

        # Fit power-law model
        self.fit_power_law()

    def apply(self, data):
        """
        Filters the catalog to select galaxies within the given redshift and stellar mass range.
        Returns a boolean mask.
        """
        return (data['z'] > self.z_min) & (data['z'] <= self.z_max) & \
               (data['SM'] > self.SM_min) & (data['SM'] <= self.SM_max)

    def measure_w_theta(self):
        """
        Computes the angular correlation function (w_theta) using CorrelationFunction class.
        """
        if self.N == 0:
            raise ValueError("No galaxies in the selected subsample. Adjust redshift/mass limits.")

        self.corr_func = CorrelationFunction(self.filtered_catalog, self.randoms, self.config)
        self.corr_func.process()

        # Extract results
        self.w_theta, self.var_w_theta, self.theta, self.rr, self.dr, self.dd = self.corr_func.calculate_w_theta()
        
         # Bootstrap errors
        self.var_w_theta_bootstrap, self.covariance_w_theta_bootstrap = self.corr_func.bootstrap_w_theta(num_bootstrap=100)

    def compute_gg(self):
        """
        Computes the matter-matter and galaxy-galaxy angular correlation functions using halomod.
        """
        if self.theta is None:
            raise ValueError("Theta values must be computed before calculating gg.")

        theta_min = np.min(self.theta) * np.pi / 180  # Convert to radians
        theta_max = np.max(self.theta) * np.pi / 180
        theta_num = len(self.theta)
        z_mean = (self.z_min + self.z_max) / 2

        self.gg = hm.integrate_corr.AngularCF(
            self.nz, self.nz,
            theta_min=theta_min,
            theta_max=theta_max,
            theta_num=theta_num,
            zmin=self.z_min,
            zmax=self.z_max,
            z= z_mean, 
            p_of_z=True
        )


        # Set HOD parameters
        self.gg.hod_params = {"M_min": 13.27, "M_1": 14.6, "alpha": 1.0}

        self.xi_g = self.gg.angular_corr_gal # galaxy-galaxy angular correlation function
        self.xi_m = self.gg.angular_corr_matter # matter-matter angular correlation function


        
        return self.gg, self.xi_g, self.xi_m

    def power_law_model(self, theta, A, gamma):

        return A * theta ** (-gamma)

    def fit_power_law(self):
        """
        Fits a power-law model to the computed w_theta values.
        Stores (A, gamma) parameters.
        """
        if self.w_theta is None or self.theta is None:
            raise ValueError("w_theta and theta must be computed before fitting a power law.")

        # Initial guesses for A and gamma
        initial_guess = [1.0, 0.8]

        print("Theta:", self.theta)
        print("w_theta:", self.w_theta)
        print("Any NaNs in theta?", np.any(np.isnan(self.theta)))
        print("Any NaNs in w_theta?", np.any(np.isnan(self.w_theta)))
        print("Any infs in w_theta?", np.any(np.isinf(self.w_theta)))


        # Fit the power-law model
        popt, _ = curve_fit(self.power_law_model, self.theta, self.w_theta, p0=initial_guess,maxfev=10000)
      
        self.power_law_params = popt  # (A, gamma)

    def get_results(self):
        """
        Returns the computed properties of the subsample:
        - Power-law parameters (A, gamma)
        - w_theta values
        - DD, DR, RR counts
        - Redshift distribution (nz)
        - gg
        """
        return {
            'N': self.N,  # Number of galaxies in the subsample
            'power_law_params': self.power_law_params,
            'w_theta': self.w_theta,
            'var_w_theta': self.var_w_theta, 
            'var_w_theta_bootstrap': self.var_w_theta_bootstrap,
            'covariance_w_theta_bootstrap': self.covariance_w_theta_bootstrap,
            'theta': self.theta,
            'dd_counts': self.dd.npairs if self.dd else None,
            'dr_counts': self.dr.npairs if self.dr else None,
            'rr_counts': self.rr.npairs if self.rr else None,
            'nz': self.nz,  
            'gg': self.gg,
            'xi_m':self.xi_m, # Matter-matter correlation
            'xi_g': self.xi_g # Galaxy-galaxy correlation
            
        }
