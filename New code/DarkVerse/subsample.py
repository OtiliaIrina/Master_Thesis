
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

        # Mask for selecting galaxies within the subsample
        self.mask = self.apply(catalog)
        self.filtered_catalog = catalog[self.mask]  # Subsample of galaxies
        
        # Initialize attributes for correlation function
        self.corr_func = None
        self.power_law_params = None  # (A, gamma) from power-law fit

        # Measure correlation function
        if w_theta is None or theta is None:
            self.measure_w_theta()
        else:
            self.w_theta = w_theta
            self.theta = theta

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
        if len(self.filtered_catalog) == 0:
            raise ValueError("No galaxies in the selected subsample. Adjust redshift/mass limits.")

        self.corr_func = CorrelationFunction(self.filtered_catalog, self.randoms, self.config)
        self.corr_func.process()

        # Extract results
        self.w_theta, self.var_w_theta, self.theta, self.rr, self.dr, self.dd = self.corr_func.calculate_w_theta()


    def power_law_model(self, theta, A, gamma):
        """
        Power-law model: w_theta = A * theta^(-gamma)
        """
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

        # Fit the power-law model
        popt, _ = curve_fit(self.power_law_model, self.theta, self.w_theta, p0=initial_guess)

        self.power_law_params = popt  # (A, gamma)

    def get_results(self):
        """
        Returns the computed properties of the subsample:
        - Power-law parameters (A, gamma)
        - w_theta values
        - DD, DR, RR counts
        """
        return {
            'power_law_params': self.power_law_params,
            'w_theta': self.w_theta,
            'theta': self.theta,
            'dd_counts': self.dd.npairs if self.dd else None,
            'dr_counts': self.dr.npairs if self.dr else None,
            'rr_counts': self.rr.npairs if self.rr else None
        }
