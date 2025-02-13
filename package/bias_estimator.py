import numpy as np
from scipy.optimize import curve_fit

class BiasEstimator:
    
    def __init__(self, gg, w_theta_list, w_theta_error_list, sum_rr_list, names, mask):
        self.gg = gg
        self.w_theta_list = w_theta_list
        self.w_theta_error_list = w_theta_error_list
        self.sum_rr_list = sum_rr_list
        self.names = names
        self.mask = mask
        self.IC_list = self.calculate_IC_values()
    
    def integral_constraint(self, wg_no_ic, rr):
        return np.sum(wg_no_ic * rr) / np.sum(rr)
    
    def calculate_IC_values(self):
        return [self.integral_constraint(self.gg.angular_corr_gal, rr) for rr in self.sum_rr_list]
    
    
    def model(self, wdm, b, IC_value):
        wg_no_ic = wdm * b**2  # Initial model without IC correction
        wg_model = wg_no_ic - IC_value  # Apply correction
        return wg_model[self.mask]  # Get only the 2-halo term
    
    def estimate_bias(self):
        b = []  # Bias values
        be = []  # Bias errors
        
        for i in range(len(self.names)):
            IC_value = self.IC_list[i]
            
            matter_corr = np.array(self.gg.angular_corr_matter).flatten()
            galaxy_corr = np.array(self.w_theta_list[i]).flatten()
            errors = np.array(self.w_theta_error_list[i]).flatten()
            
            min_len = min(len(matter_corr), len(galaxy_corr), len(errors))
            matter_corr = matter_corr[:min_len]
            galaxy_corr = galaxy_corr[:min_len]
            errors = errors[:min_len]
            
            popt, pcov = curve_fit(lambda wdm, b: self.model(wdm, b, IC_value),
                                   matter_corr, galaxy_corr[self.mask],
                                   sigma=errors[self.mask], absolute_sigma=True, p0=[1.0])
            
            bias = popt[0]
            bias_error = np.sqrt(np.diag(pcov))[0]
            
            b.append(bias)
            be.append(bias_error)
            print(f"{self.names[i]}: Bias = {bias:.3f} ± {bias_error:.3f}")
        
        return b, be
