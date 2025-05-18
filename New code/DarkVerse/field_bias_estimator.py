import numpy as np
from scipy.optimize import curve_fit

class FieldBiasEstimator:
    def __init__(self, xi_m, mask):
        """
        Initialize with matter correlation and mask that's common to all fields
        
        Parameters:
        - xi_m: np.array, matter-matter angular correlation function
        - mask: np.array, boolean mask for selecting the 2-halo term
        """
        self.xi_m = xi_m
        self.mask = mask
        self.fields = []
        
    def add_field(self, field_name, xi_g, w_theta, w_theta_error, sum_rr):
        """
        Add data for a single field
        
        Parameters:
        - field_name: str, identifier for the field
        - xi_g: np.array, galaxy-galaxy correlation for this field
        - w_theta: np.array, observed correlation function
        - w_theta_error: np.array, errors for observed correlation
        - sum_rr: np.array, RR pair counts for this field
        """
        self.fields.append({
            'name': field_name,
            'xi_g': xi_g, 
            'w_theta': w_theta,
            'w_error': w_theta_error,
            'sum_rr': sum_rr,
            'IC': None,
            'bias': None,
            'bias_error': None
        })
    
    def _calculate_IC(self, wg_no_ic, sum_rr):
        """Calculate integral constraint for given correlation and RR counts"""
        return np.sum(wg_no_ic * sum_rr) / np.sum(sum_rr)
    
    def _model(self, wdm, b, IC):
        """Model function with current IC"""
        return (wdm * b**2 - IC)[self.mask]
    
    def fit_individual_fields(self):
        """Fit bias for each field separately"""
        results = []
        for field in self.fields:
            try:
                # Initial IC estimate using xi_g
                field['IC'] = self._calculate_IC(field['xi_g'], field['sum_rr'])
                
                # Perform fit
                popt, pcov = curve_fit(
                    lambda wdm, b: self._model(wdm, b, field['IC']),
                    self.xi_m,
                    field['w_theta'][self.mask], 
                    sigma=field['w_error'][self.mask],
                    p0=[1.0],
                    maxfev=10000
                )
                
                field['bias'] = popt[0]
                field['bias_error'] = np.sqrt(np.diag(pcov))[0]
                
                # Update IC with best-fit model
                wg_no_ic_fit = self.xi_m * field['bias']**2
                field['IC'] = self._calculate_IC(wg_no_ic_fit, field['sum_rr'])
                
                results.append((field['name'], field['bias'], field['bias_error'], field['IC']))
                
            except RuntimeError:
                results.append((field['name'], np.nan, np.nan, np.nan))
        
        return results
    
    def fit_combined(self):
        """Fit combined data from all fields to get global bias"""
        if not self.fields:
            raise ValueError("No fields added for fitting")
            
        # Stack all field data
        all_w_theta = np.array([f['w_theta'][self.mask] for f in self.fields])
        all_errors = np.array([f['w_error'][self.mask] for f in self.fields])
        weights = 1.0 / (all_errors**2)
        
        # Weighted average of the observations
        combined_w = np.sum(all_w_theta * weights, axis=0) / np.sum(weights, axis=0)
        combined_err = 1.0 / np.sqrt(np.sum(1.0 / all_errors**2, axis=0))
        
        # Use mean of individual ICs as starting point
        initial_IC = np.mean([f['IC'] for f in self.fields if f['IC'] is not None])
        
        try:
            popt, pcov = curve_fit(
                lambda wdm, b: self._model(wdm, b, initial_IC),
                self.xi_m,
                combined_w,
                sigma=combined_err,
                p0=[1.0],
                maxfev=10000
            )
            
            global_bias = popt[0]
            global_error = np.sqrt(np.diag(pcov))[0]
            
            return global_bias, global_error
            
        except RuntimeError:
            return np.nan, np.nan

    def get_field_models(self, global_bias):
        models = {}
        for field in self.fields:
            if field['IC'] is not None:
                models[field['name']] = self._model(self.xi_m, global_bias, field['IC'])
            else:
                models[field['name']] = None
        return models
             