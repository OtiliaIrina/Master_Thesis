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

# Usage example:
high_sm_biases = []
for subsample in subsamples:
    # Create a BiasCalculator object
    calculator = BiasCalculator(w_model, gg)

    # Calculate bias and error for the subsample
    bias, bias_error, popt, pcov = calculator.calculate_bias(subsample.info['w'])

    # Print and store results
    print(f"Subsample {i+1} SM_mean: {subsample.info['SM_mean']:.3f}")
    print(f"Bias: {bias:.3f} +/- {bias_error:.3f}")

    high_sm_biases.append([subsample.info['SM_mean'], bias, bias_error])

# ... rest of your code for saving and loading data using high_sm_biases