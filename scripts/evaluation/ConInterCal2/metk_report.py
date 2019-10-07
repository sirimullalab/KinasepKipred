#!/usr/bin/env python

from __future__ import print_function
import sys
import pandas as pd
from metk_util import rmse, mean_absolute_error, pearson_confidence,spearman_confidence,max_possible_correlation, ci_confidence, rmse_confidence
from scipy.stats import pearsonr, kendalltau, spearmanr
from evaluation_metrics import ci
def metk_report(df_pki_cal):
    """
    Generate a report
    :param df_pki_cal: input dataframe, with actual and predicted results
    :param outfile: output file for the report
    :return: the report as a list of strings
    """


    N = df_pki_cal.shape[0]
    pred = df_pki_cal['Pred']
    expr = df_pki_cal['Exp']
#    rms_val = rmse(pred, expr)
#    mae_val = mean_absolute_error(pred, expr)
    ci_score = ci(expr, pred)
    rms = rmse(pred, expr)
    pearson_r, pearson_p = pearsonr(pred, expr)
    spearman_r, spearman_p = spearmanr(pred, expr)
#    kendall_t, kendall_p = kendalltau(pred, expr)
    pearson_vals = [x for x in [pearson_r] + list(pearson_confidence(pearson_r, N))]
    spearman_vals = [x for x in [spearman_r] + list(spearman_confidence(spearman_r, N))]
    rmse_vals = [x for x in [rms] + list(rmse_confidence(rms, N))]
    ci_vals = [x for x in [ci_score] + list(ci_confidence(ci_score, N))]

    report = []
    report.append("N = %d" % N)
#    report.append("RMSE = %.2f kcal/mol" % rms_val)
#    report.append("MAE  = %.2f kcal/mol" % mae_val)
    report.append("Pearson = %0.3f  95%%ConInterval = %.3f %.3f" % tuple(pearson_vals))
    report.append("Spearman = %0.3f 95%%ConInterval = %.3f %.3f" % tuple(spearman_vals))
    report.append("rmse = %0.3f 95%%ConInterval = %.3f %.3f" % tuple(rmse_vals))
#    report.append("Kendall tau = %0.2f" % kendall_t)
    report.append("ci = %0.3f  95%%ConInterval = %.3f %.3f" % tuple(ci_vals))

    return report

def main():
    df_pki_cal = pd.read_csv(sys.argv[1])
    metk_report(df_pki_cal)

if __name__ == "__main__":
    main()
