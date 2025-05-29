from .correlation_function import CorrelationFunction
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import least_squares

import halomod as hm



class Selection:
    def __init__(self, catalog, randoms, z_min, z_max, SM_min, SM_max, config, w_theta=None, theta=None):
        """
        Initialize a subsample of galaxies with given redshift and stellar mass limits.

        Parameters:
        - catalog: Galaxy catalog data
        - randoms: Random catalog for correlation function
        - z_min, z_max: Redshift range limits
        - SM_min, SM_max: Stellar mass range limits
        - config: Configuration for correlation function computation
        - w_theta, theta: Optional pre-computed correlation function
        """
        self.catalog = catalog
        self.randoms = randoms
        self.z_min = z_min
        self.z_max = z_max
        self.SM_min = SM_min
        self.SM_max = SM_max
        self.config = config
        self.info = {}  # Dictionary to store additional information

        # Compute the redshift distribution
        self.nz = hm.integrate_corr.flat_z_dist(self.z_min, self.z_max)

        # Apply selection to create subsample
        self.mask = self.apply(self.catalog)
        self.filtered_catalog = self.catalog[self.mask]
        self.N = len(self.filtered_catalog)

        # Initialize correlation function attributes
        self.corr_func = None
        self.power_law_params = None
        self.gg = None
        self.xi_g = None
        self.xi_m = None
        self.hod_params = None

        # Measure or set correlation function
        if w_theta is None or theta is None:
            self.measure_w_theta()
        else:
            self.w_theta = w_theta
            self.theta = theta

        # Initialize with default HOD (no fitting yet)
        self._init_halomod()

        # Fit power-law model
        self.fit_power_law()

    def apply(self, data):
        """Filter catalog to select galaxies within redshift and mass limits."""
        return (
            (data['z'] > self.z_min) & (data['z'] <= self.z_max) &
            (data['SM'] > self.SM_min) & (data['SM'] <= self.SM_max)
        )

    def measure_w_theta(self):
        """Compute the angular correlation function from the data."""
        if self.N == 0:
            raise ValueError("No galaxies in the selected subsample.")

        self.corr_func = CorrelationFunction(self.filtered_catalog, self.randoms, self.config)
        self.corr_func.process()

        # Extract correlation function results
        (
            self.w_theta,
            self.var_w_theta,
            self.theta,
            self.rr,
            self.dr,
            self.dd
        ) = self.corr_func.calculate_w_theta()

        # Compute bootstrap errors
        (
            self.var_w_theta_bootstrap,
            self.covariance_w_theta_bootstrap
        ) = self.corr_func.bootstrap_w_theta(num_bootstrap=100)

    def power_law_model(self, theta, A, gamma):
        """Power-law model for correlation function."""
        return A * theta ** (-gamma)

    def fit_power_law(self):
        """Fit power-law model to the measured correlation function."""
        if self.w_theta is None or self.theta is None:
            raise ValueError("w_theta and theta must be computed first.")

        initial_guess = [1.0, 0.8]
        popt, _ = curve_fit(
            self.power_law_model,
            self.theta,
            self.w_theta,
            p0=initial_guess,
            maxfev=10000
        )
        self.power_law_params = popt  # (A, gamma)

    def _init_halomod(self):
        """Initialize halo model with default parameters."""
        if self.theta is None:
            raise ValueError("Theta values must be computed first.")

        theta_min = np.min(self.theta) * np.pi / 180
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
            z=z_mean,
            p_of_z=True
        )

        self.gg.hod_params =  {"M_min": 12.0, "M_1": 13.0, "alpha": 1.0, "central":True} #{"M_min": 12.5, "M_1": 13.5, "alpha": 1.0}
        

        self.xi_g = self.gg.angular_corr_gal
        self.xi_m = self.gg.angular_corr_matter

        return self.gg, self.xi_g, self.xi_m

    def hod_model(self, M_min, M_1, alpha):
        """
        Compute galaxy correlation function for given HOD parameters.

        Parameters:
        - M_min: log10 of minimum halo mass [log10(Msun/h)]
        - M_1: log10 of characteristic halo mass [log10(Msun/h)]
        - alpha: power law slope for satellites

        Returns:
        - xi_g: model galaxy correlation function
        """
       

        self.gg.hod_params = {
            "M_min": M_min,
            "M_1": M_1,
            "alpha": alpha
        }

        self.gg.update()  # Try this ?

        """ Includ the integral constrain bias in the hod model that we fit"""
        self.IC = np.sum(self.gg.angular_corr_gal * self.rr.npairs) / np.sum(self.rr.npairs)

        return self.gg.angular_corr_gal - self.IC

    def fit_hod(self, theta=None, p0=None, bounds=None):
        """
        Fit 3-parameter HOD model to observed correlation function using curve_fit.

        Parameters:
        - theta: Angular scales to fit 
        - p0: Initial guess [logM_min, logM_1, alpha]
        - bounds

        Returns:
        - popt: Parameters [logM_min, logM_1, alpha]
        - pcov: Parameter covariance matrix
        """

        
        #if p0 is None:
           # p0 = [12.5, 13.5, 1.0]
        #if bounds is None:
           # bounds = ([11.0, 12.5, 0.3], [14.5, 15.5, 2.0]) # bounds = ([12.0, 13.0, 0.5], [13.0, 14.0, 1.5])
        

        """We can improve the sensitivity of the M_min parameter! 
        Curve fit takes a tiny step to begin exploring, so we need to make its steps more meaningful.
        We can do that by passing it a much smaller value, and then mutliplying it when using it. E.g."""
        
        if p0 is None:
            p0 = [12.5*1e-7, 13.5, 1.0]
        if bounds is None:
            bounds = ([11.0*1e-7, 12.5, 0.3], [14.5*1e-7, 15.5, 2.0]) 


        theta = self.theta if theta is None else theta   
        w = self.w_theta if theta is None else np.interp(theta, self.theta, self.w_theta)
        err = np.sqrt(self.var_w_theta_bootstrap)

        # Fit only the 2-halo term (only 0.1 <= theta <= 0.4 degrees)
        mask_theta = (theta >= 0.1) & (theta <= 0.4)
        theta_fit = theta[mask_theta]
        w_fit = w[mask_theta]
        err_fit = err[mask_theta]

        
        """And then multiply that big factor, 1e7, when passing it to the HOD package,
        self.hod_model(logM_min*1e7, logM_1, alpha)[mask_theta]"""

        def hod_wrapper(theta, logM_min, logM_1, alpha):
            return self.hod_model(logM_min*1e7, logM_1, alpha)[mask_theta] #self.hod_model(logM_min, logM_1, alpha)[mask_theta]

        self.hod_params, pcov = curve_fit(
            hod_wrapper,
            theta_fit, w_fit,
            p0=p0,
            sigma=err_fit,
            bounds=bounds,
            maxfev=10000
        )

        return self.hod_params, pcov


        # trying to "unstuck" Mmin aka improve the sensitivity
        """
        def residuals(params):
            logM_min, logM_1, alpha = params
            model = self.hod_model(logM_min, logM_1, alpha)
            return (w - model) / err

        result = least_squares(
            residuals,
            p0,
            bounds=bounds,
            max_nfev=10000,
            jac='2-point',
            verbose=0
        )


        self.hod_params = result.x

        pdict = {
            "residuals": result.fun,
            "jacobian": result.jac,
            "cost": result.cost,
            "success": result.success,
            "message": result.message,
            "nfev": result.nfev
        }

        return self.hod_params, None, pdict"""

    def fit_shmr_fixed_alpha(self, theta=None, p0=None, bounds=None):
        """
        Fit HOD model with fixed alpha=1.0 to get better constraints on M_min,
        and compute the stellar-to-halo mass ratio (SHMR = M*/M_halo).
    
        Parameters:
        - theta: Angular scales to fit 
        - p0: Initial guess [logM_min, logM_1]
        - bounds: Bounds for [logM_min, logM_1]
    
        Returns:
        - params: Fitted parameters [logM_min, logM_1, fixed alpha]
        - shmr: Estimated SHMR = 10^mean stellar mass / 10^M_min
        """
        fixed_alpha = 1.0
    
        if p0 is None:
            p0 = [12.5 * 1e-7, 13.5]
        if bounds is None:
            bounds = ([11.0 * 1e-7, 12.5], [14.5 * 1e-7, 15.5])
    
        theta = self.theta if theta is None else theta
        w = self.w_theta if theta is None else np.interp(theta, self.theta, self.w_theta)
        err = np.sqrt(self.var_w_theta_bootstrap)
    
        # Use only large-scale (2-halo) part of the correlation function
        mask_theta = (theta >= 0.1) & (theta <= 0.4)
        theta_fit = theta[mask_theta]
        w_fit = w[mask_theta]
        err_fit = err[mask_theta]
    
        def hod_wrapper_fixed_alpha(theta_vals, logM_min, logM_1):
            return self.hod_model(logM_min * 1e7, logM_1, fixed_alpha)[mask_theta]
    
        popt, _ = curve_fit(
            hod_wrapper_fixed_alpha,
            theta_fit, w_fit,
            p0=p0,
            sigma=err_fit,
            bounds=bounds,
            maxfev=10000
        )
    
        # Store fitted params
        logM_min_fit, logM_1_fit = popt
        self.hod_params = [logM_min_fit * 1e7, logM_1_fit, fixed_alpha]
    
        # Compute SHMR = M_star / M_halo
        mean_stellar_mass_log10 = 0.5 * (self.SM_min + self.SM_max)
        shmr = 10**mean_stellar_mass_log10 / 10**(logM_min_fit * 1e7)
    
        return self.hod_params, shmr
    

    def get_results(self):
        """
        Return dictionary with all computed results.

        Includes:
        - Basic sample properties
        - Correlation function measurements
        - Power law fit parameters
        - HOD fit parameters
        - Model correlation functions
        """
        return {
            'N': self.N,
            'z_range': (self.z_min, self.z_max),
            'SM_range': (self.SM_min, self.SM_max),

            # Correlation function measurements
            'w_theta': self.w_theta,
            'theta': self.theta,
            'var_w_theta': self.var_w_theta,
            'var_w_theta_bootstrap': self.var_w_theta_bootstrap,
            'covariance_w_theta_bootstrap': self.covariance_w_theta_bootstrap,

            # Pair counts
            'dd_counts': self.dd.npairs if hasattr(self, 'dd') and self.dd else None,
            'dr_counts': self.dr.npairs if hasattr(self, 'dr') and self.dr else None,
            'rr_counts': self.rr.npairs if hasattr(self, 'rr') and self.rr else None,

            # Power law fit
            'power_law_params': self.power_law_params,

            # HOD model
            'hod_params': self.hod_params,
            'shmr': ( 10**(0.5 * (self.SM_min + self.SM_max)) / 10**self.hod_params[0] if self.hod_params is not None else None),


            # Model correlations
            'nz': self.nz,
            'xi_m': self.xi_m,
            'xi_g': self.xi_g,
            'gg_model': self.gg
        }
