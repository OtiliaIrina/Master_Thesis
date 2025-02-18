import treecorr
import numpy as np
from scipy.optimize import curve_fit
import halomod as hm

class Subsample:
    def __init__(self, catalog, z_min, z_max, SM_min, SM_max, w_theta=None, theta=None):
        self.catalog = catalog
        self.z_min = z_min
        self.z_max = z_max
        self.SM_min = SM_min
        self.SM_max = SM_max
        self.w_theta = w_theta  # Observed angular correlation function
        self.theta = theta  # Corresponding angle values (in degrees)
        
        # Mask for selecting galaxies within the subsample
        self.mask = self.apply(catalog)
        
        # Extract relevant data
        self.filtered_catalog = catalog[self.mask]  # Subsample of galaxies
        
        # Initialize attributes for correlation function
        self.dd = None
        self.dr = None
        self.rr = None
        self.power_law_params = None  # (A, gamma) from power-law fit

        # Measure correlation function if data is available
        if self.w_theta is None or self.theta is None:
            self.measure_w_theta()
        
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
        Computes the angular correlation function (w_theta) using TreeCorr.
        Stores DD, DR, RR counts.
        """
        config = {
            'min_sep': 0.01,  # Minimum separation in degrees
            'max_sep': 1.0,   # Maximum separation in degrees
            'nbins': 10,      # Number of bins
            'sep_units': 'degrees'
        }

        # Generate TreeCorr catalogs
        data_cat = treecorr.Catalog(ra=self.filtered_catalog['ra'], dec=self.filtered_catalog['dec'],
                                    ra_units='degrees', dec_units='degrees', npatch=50)
        
        # Generate random catalog for comparison
        random_cat = treecorr.Catalog(ra=np.random.uniform(0, 360, len(self.filtered_catalog)),
                                      dec=np.random.uniform(-90, 90, len(self.filtered_catalog)),
                                      ra_units='degrees', dec_units='degrees', npatch=50)

        # Create correlation function estimators
        self.dd = treecorr.NNCorrelation(config)
        self.dr = treecorr.NNCorrelation(config)
        self.rr = treecorr.NNCorrelation(config)

        # Process the data
        self.dd.process(data_cat)
        self.dr.process(data_cat, random_cat)
        self.rr.process(random_cat)

        # Compute w(theta)
        self.theta = np.exp(self.dd.meanlogr)
        self.w_theta, _ = self.dd.calculateXi(rr=self.rr, dr=self.dr)

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
