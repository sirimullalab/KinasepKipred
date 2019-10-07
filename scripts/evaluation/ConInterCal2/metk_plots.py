#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from metk_util import check_dataframe
import math

def add_kcal_error(df, bins=None):
    """
    Add columns to a dataframe showing absolute and binned err
    :param df: input dataframe
    :param bins: bins to use (currently "<1 kcal", "1-2 kcal", ">2 kcal")
    :return: None
    """
    if bins is None:
        bins = [1, 2]
    pt_color = ['green', 'yellow', 'red']
    df['Error'] = np.abs(df['Exp'] - df['Pred'])
    df['Error_Bin'] = [pt_color[x] for x in np.digitize(df['Error'], bins)]

def kcal_plot(df, ax, axis_range=None):
    """
    Draw a scatterplot of experimental vs predicted Ki or IC50
    :param df: input dataframe
    :param ax: matplotlib axis
    :param axis_range: range for axes [minX, maxY, minY, maxY]
    :return: None
    """
    if axis_range is None:
        axis_range = [-12, -6, -12, -6]
    pt_color = ['green', 'yellow', 'red']
    df['Error'] = np.abs(df['Exp'] - df['Pred'])
    df['Error_Bin'] = [pt_color[x] for x in np.digitize(df['Error'], [1, 2])]
    ax.axis(axis_range)
    ax.set_xlabel(r'Actual pKi')
    ax.set_ylabel(r'Predicted pKi')
    ax.scatter(df['Exp'], df['Pred'], s=100, c=df['Error_Bin'], alpha=0.5, edgecolors="black")
    # y = x
    ax.plot([0, -100], [0, -100], linewidth=2, color='black')
    # y = x+1
    ax.plot([-1, -100], [0, -100], linewidth=1, color="blue", linestyle='--')
    # y = x-1
    ax.plot([1, -100], [0, -100], linewidth=1, color="blue", linestyle='--')
    # y = x+2
    ax.plot([-2, -100], [0, -100], linewidth=1, color="black")
    # y = x-2
    ax.plot([2, -100], [0, -100], linewidth=1, color="black")


def draw_plots(df_kcal,pdf_file_name):
    """
    Draw scatter plots and histograms showing agreement between experimental and predicted activity
    :param df_kcal: input dataframe, data is in kcal/mol
    :param pdf_file_name: output file for plot
    :param units: units to use for the plots (currently uM or nM)
    :return:
    """
    add_kcal_error(df_kcal)
    f_kcal, ax_kcal = plt.subplots(1, figsize=(7, 7))
    ax_kcal.set_title("N = %d" % df_kcal.shape[0])
    
    minx = int( min(df_kcal["Exp"] ) - 1 )
    maxx = int( max(df_kcal["Exp"] ) + 1 )
    miny = int( min(df_kcal["Pred"]) - 1 )
    maxy = int( max(df_kcal["Pred"]) + 1 )
    
    kcal_plot(df_kcal, ax_kcal, axis_range=[minx, maxx, miny, maxy])
    
    pdf_pages = PdfPages(pdf_file_name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.tight_layout()
    pdf_pages.savefig(f_kcal.get_figure())

    pdf_pages.close()

def generate_pdf(figure_list, pdf_file_name):
    pdf_pages = PdfPages(pdf_file_name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.tight_layout()
    for ax in figure_list:
        pdf_pages.savefig(ax.get_figure())
    pdf_pages.close()


def main():
    pdf_file_name = "/Users/gvin/pki_paper_work/Kinase_pKI/metk/modelevaltoolkit/collect/results_filename4.pdf" 
    df_kcal = pd.read_csv(sys.argv[1])
    check_dataframe(df_kcal)
    draw_plots(df_kcal, pdf_file_name)


if __name__ == "__main__":
    main()
