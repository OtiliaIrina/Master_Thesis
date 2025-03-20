import numpy as np
from scipy.optimize import curve_fit


class BiasEstimator:
    
    def __init__(self, xi_g, xi_m, w_theta_list, w_theta_error_list, sum_rr_list, mask):
        """
        Initializes the BiasEstimator class.

        Parameters:
        - xi_g: np.array, galaxy-galaxy angular correlation function
        - xi_m: np.array, matter-matter angular correlation function
        - w_theta_list: list, observed angular correlation functions
        - w_theta_error_list: list, corresponding errors
        - sum_rr_list: list, RR pair counts
        - mask: np.array, boolean mask for selecting the 2-halo term
        """
        self.xi_g = xi_g  # Galaxy-galaxy angular correlation function      
        self.xi_m = xi_m  # Matter-matter angular correlation function
        self.w_theta_list = w_theta_list
        self.w_theta_error_list = w_theta_error_list
        self.sum_rr_list = sum_rr_list
        self.mask = mask
        self.IC_list = self.calculate_IC_values()  # Compute Integral Constraint (IC)

    def integral_constraint(self, wg_no_ic, rr):
        """
        Computes the integral constraint correction.

        Parameters:
        - wg_no_ic: np.array, galaxy correlation without IC correction
        - rr: np.array, RR pair counts

        Returns:
        - IC value
        """
        return np.sum(wg_no_ic * rr) / np.sum(rr)
    
    def calculate_IC_values(self):
        """
        Computes the integral constraint correction for each subsample.

        Returns:
        - List of IC values for each subsample
        """
        return [self.integral_constraint(self.xi_g, rr) for rr in self.sum_rr_list]
    
    def model(self, wdm, b, IC_value):
        """
        Model function for bias estimation.

        Parameters:
        - wdm: np.array, matter-matter correlation function
        - b: float, bias parameter
        - IC_value: float, Integral Constraint (IC)

        Returns:
        - Modeled galaxy correlation function
        """
        wg_no_ic = wdm * b**2  # Initial model without IC correction
        wg_model = wg_no_ic - IC_value  # Apply correction
        return wg_model[self.mask]  # Get only the selected range
    
    def estimate_bias(self):
        """
        Estimates the bias (b) by fitting the matter correlation function to the observed galaxy correlation.

        Returns:
        - b: List of bias values for each sample
        - be: List of bias errors
        """
        b = []  # Bias values
        be = []  # Bias errors
        
        for i in range(len(self.w_theta_list)):
            IC_value = self.IC_list[i]  # Integral Constraint for this sample
            
            # Flatten arrays for safety
            matter_corr = np.array(self.xi_m).flatten()
            galaxy_corr = np.array(self.w_theta_list[i]).flatten()
            errors = np.array(self.w_theta_error_list[i]).flatten()
            
            # Ensure arrays have the same length
            min_len = min(len(matter_corr), len(galaxy_corr), len(errors))
            matter_corr = matter_corr[:min_len]
            galaxy_corr = galaxy_corr[:min_len]
            errors = errors[:min_len]
            
            try:
                popt, pcov = curve_fit(
                    lambda wdm, b: self.model(wdm, b, IC_value),
                    matter_corr, 
                    galaxy_corr[self.mask],  # Apply mask
                    sigma=errors[self.mask], 
                    absolute_sigma=True, 
                    p0=[1.0],  # Initial guess for bias
                    maxfev=10000)  
                    
                bias = popt[0]
                bias_error = np.sqrt(np.diag(pcov))[0]
                
                b.append(bias)
                be.append(bias_error)
            
            except RuntimeError:
                b.append(np.nan)
                be.append(np.nan)

        return b, be
