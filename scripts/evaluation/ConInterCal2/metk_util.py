#!/usr/bin/env python

from __future__ import print_function
import sys
import math
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import norm
from scipy.stats import spearmanr
from evaluation_metrics import AUC, average_AUC

def rmse(pred_array, ref_array):
    """
    Calculate root mean squared (rms) error
    :param pred_array: the predicted values
    :param ref_array: the reference values
    :return: the rms error
    """
    return np.sqrt(np.mean((pred_array - ref_array) ** 2))


def mean_absolute_error(pred_array, ref_array):
    """
    Calculate mean absolute error
    :param pred_array: the predicted values
    :param ref_array: the reference values
    :return: the mean absolute error
    """
    return np.mean(np.abs(pred_array - ref_array))


def get_unit_multiplier(units):
    """
    Function so that I only have to put the unit dictionary in one place
    :param: units: units
    :return: unit dictionary
    """
    multiplier_dict = {"M": 1, "mM": 1e-3, "uM": 1e-6, "nM": 1e-9}
    try:
        multiplier = multiplier_dict[units]
        return multiplier
    except KeyError:
        print("Error:", units, "is not supported in ki_to_kcal")
        sys.exit(0)


def pearson_confidence(r, num, interval=0.95):
    """
    Calculate upper and lower 95% CI for a Pearson r
    Inspired by https://stats.stackexchange.com/questions/18887
    :param r: Pearson's R
    :param num: number of data points
    :param interval: confidence interval (0-1.0)
    :return: lower bound, upper bound
    """
    stderr = 1.0 / math.sqrt(num - 3)
    z_score = norm.ppf(interval)
    delta = z_score * stderr
    lower = math.tanh(math.atanh(r) - delta)
    upper = math.tanh(math.atanh(r) + delta)
    return lower, upper

def spearman_confidence(rho, num, interval=0.95):
    """ 
    Calculate upper and lower 95% CI for a spearman (not R**2)
    Inspired by https://stats.stackexchange.com/questions/18887
    :param r: spearman's rho
    :param num: number of data points
    :param interval: confidence interval (0-1.0)
    :return: lower bound, upper bound
    """
    stderr = 1.0 / math.sqrt(num - 3)
    z_score = norm.ppf(interval)
    delta = z_score * stderr
    lower = math.tanh(math.atanh(rho) - delta)
    upper = math.tanh(math.atanh(rho) + delta)
    return lower, upper
def ci_confidence(ci_score, num, interval=0.95):
    """ 
    Calculate upper and lower 95% CI for a Concordance Interval)
    Inspired by https://stats.stackexchange.com/questions/18887
    :param r: spearman's rho
    :param num: number of data points
    :param interval: confidence interval (0-1.0)
    :return: lower bound, upper bound
    """
    stderr = 1.0 / math.sqrt(num - 3)
    z_score = norm.ppf(interval)
    delta = z_score * stderr
    lower = math.tanh(math.atanh(ci_score) - delta)
    upper = math.tanh(math.atanh(ci_score) + delta)
    return lower, upper
def rmse_confidence(rmse_value, num, interval=0.95):
    """ 
    Calculate upper and lower 95% CI for root mean squre error
    Inspired by https://stats.stackexchange.com/questions/18887
    :param r: spearman's rho
    :param num: number of data points
    :param interval: confidence interval (0-1.0)
    :return: lower bound, upper bound
    """
    stderr = 1.0 / math.sqrt(num - 3)
    z_score = norm.ppf(interval)
    delta = z_score * stderr
    lower = math.tanh(math.atanh(rmse_value) - delta)
    upper = math.tanh(math.atanh(rmse_value) + delta)
    return lower, upper

def AUC_confidence(auc_value, num, interval=0.95):
    """ 
    Calculate upper and lower 95% CI for area under the roc curve
    Inspired by https://stats.stackexchange.com/questions/18887
    :param r: spearman's rho
    :param num: number of data points
    :param interval: confidence interval (0-1.0)
    :return: lower bound, upper bound
    """
    stderr = 1.0 / math.sqrt(num - 3)
    z_score = norm.ppf(interval)
    delta = z_score * stderr
    lower = math.tanh(math.atanh(auc_value) - delta)
    upper = math.tanh(math.atanh(auc_value) + delta)
    return lower, upper
def max_possible_correlation(vals, error=1 / 3.0, method=pearsonr, cycles=1000):
    """
    Calculate the maximum possible correlation given a particular experimental error
    Based on Brown, Muchmore, Hajduk http://www.sciencedirect.com/science/article/pii/S1359644609000403
    :param vals: experimental values (should be on a log scale)
    :param error: experimental error
    :param method: method for calculating the correlation, must take 2 lists and return correlation and p_value
    :param cycles: number of random cycles
    :return: maximum possible correlation
    """
    cor_list = []
    for i in range(0, cycles):
        noisy_vals = []
        for val in vals:
            noisy_vals.append(val + np.random.normal(0, error))
        cor_list.append(method(vals, noisy_vals)[0])
    return np.mean(cor_list)



def check_dataframe(df):
    """
    Check a dataframe to ensure that we have columns called "Pred" and "Exp"
    :param df: input dataframe
    :return: no return value, exits on error
    """
    cols = df.columns
    if "Pred" not in cols:
        print('Input Error: Your input file does not have a column named "Pred"', file=sys.stderr)
        sys.exit(0)
    if "Exp" not in cols:
        print('Input Error: Your input file does not have a column named "Exp"', file=sys.stderr)
        sys.exit(0)
