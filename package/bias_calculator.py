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


    def __init__(self, w_model, gg):

        self.w_model = w_model
        self.gg = gg

    def calculate_bias(self, w_data):

        popt, pcov = curve_fit(self.w_model, self.gg.angular_corr_matter, w_data[:len(self.gg.angular_corr_matter)], p0=[2.0])
        bias = popt[0]
        bias_error = np.sqrt(pcov[0, 0])

        return bias, bias_error, popt, pcov
