import astropy.io.fits as fits
import numpy as np
from scipy.integrate import quad
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table,join
import os
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM



class BiasCalculator:
    """
    Calculates the bias parameter and its error from a w(theta) function and galaxy-galaxy correlation data.
    """

    def __init__(self, w_model, gg):
        """
        Initializes the BiasCalculator object.

        Args:
            w_model (function): The function representing the w(theta) model.
            gg (object): An object containing the galaxy-galaxy correlation data (e.g., angular_corr_matter).
        """
        self.w_model = w_model
        self.gg = gg

    def calculate_bias(self, w_data):
        """
        Calculates the bias parameter and its error from the provided w_data.

        Args:
            w_data (array): An array containing the w(theta) data.

        Returns:
            tuple: A tuple containing the best-fit bias, its error, and the fitting parameters.
        """

        popt, pcov = curve_fit(self.w_model, self.gg.angular_corr_matter, w_data[:len(self.gg.angular_corr_matter)], p0=[2.0])
        bias = popt[0]
        bias_error = np.sqrt(pcov[0, 0])

        return bias, bias_error, popt, pcov

