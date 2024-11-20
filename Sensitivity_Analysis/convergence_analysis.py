import numpy as np
from scipy.stats import wasserstein_distance, shapiro

def mean_variance_convergence(data):
    """
    Calculate the running mean and variance of a dataset.

    This function computes the running mean and variance for each point in the input data array.
    The running mean is calculated using the cumulative sum, and the running variance is calculated
    using the sample variance formula with Bessel's correction (ddof=1).

    Parameters:
    data (array-like): The input data for which the running mean and variance are to be calculated.

    Returns:
    tuple: A tuple containing two elements:
        - running_means (numpy.ndarray): An array of running means.
        - running_vars (list): A list of running variances.
    """
    running_means = np.cumsum(data) / np.arange(1, len(data) + 1)
    running_vars = [np.var(data[:i+1], ddof=1) for i in range(len(data))]
    return running_means, running_vars

def wasserstein_convergence(data, window_size):
    """
    Calculate the Wasserstein distance convergence over a sliding window.

    This function computes the Wasserstein distance between a reference distribution
    and a current distribution over a sliding window of specified size. The reference
    distribution grows incrementally with each step, while the current distribution
    is a fixed-size window that slides over the data.

    Parameters:
    data (list or array-like): The input data over which the Wasserstein distance is calculated.
    window_size (int): The size of the sliding window for the current distribution.

    Returns:
    list: A list of Wasserstein distances for each window position.
    """
    distances = []
    for i in range(window_size, len(data) + 1):
        reference_distribution = data[:i]
        current_distribution = data[i - window_size:i]
        distance = wasserstein_distance(reference_distribution, current_distribution)
        distances.append(distance)
    return distances

def test_convergence(data, method='auto', window_size=5):
    """
    Test convergence of replicates using either mean/variance stabilization or Wasserstein distance.

    Parameters:
        data (np.array): Array of model outputs from multiple replicates.
        method (str): Method to use for convergence testing ('auto', 'mean', 'wasserstein').
                      'auto' selects method based on normality test.
        window_size (int): Sliding window size for Wasserstein distance (only used if method is 'wasserstein').

    Returns:
        list: Convergence metric values at each step.
    """
    # Automatically select method based on normality test
    if method == 'auto':
        stat, p_value = shapiro(data)
        method = 'mean' if p_value > 0.05 else 'wasserstein'
        print(f"Auto-selected method: {'Mean/Variance' if method == 'mean' else 'Wasserstein Distance'} (p-value: {p_value:.3f})")

    # Use selected method to check convergence
    if method == 'mean':
        means, variances = mean_variance_convergence(data)
        # Convergence metric: Combine changes in means and variances
        convergence_metric = np.abs(np.diff(means)) + np.abs(np.diff(variances))
    elif method == 'wasserstein':
        convergence_metric = wasserstein_convergence(data, window_size)
    else:
        raise ValueError("Invalid method. Choose 'auto', 'mean', or 'wasserstein'.")

    return convergence_metric, method


def bootstrap_convergence(data, num_bootstrap=1000, statistic=np.mean):
    """
    Test convergence using bootstrapping by resampling the data and analyzing the stability of a statistic.

    Parameters:
        data (np.array): Array of model outputs from multiple replicates.
        num_bootstrap (int): Number of bootstrap resamples.
        statistic (function): Statistic to compute on each bootstrap sample (e.g., np.mean, np.median).

    Returns:
        list: Bootstrapped convergence metric values (confidence intervals for the statistic).
    """
    bootstrap_metrics = []
    for i in range(1, len(data) + 1):
        current_data = data[:i]
        bootstrap_samples = [statistic(np.random.choice(current_data, size=len(current_data), replace=True)) for _ in range(num_bootstrap)]
        metric = np.std(bootstrap_samples)  # Variability of the bootstrapped statistic
        bootstrap_metrics.append(metric)
    return bootstrap_metrics
