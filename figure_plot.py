#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import zsx_some_tools as st
from scipy import stats
import os
import sys
import argparse
import warnings
warnings.filterwarnings("ignore")


def plot_system_pool(p_data, lw=2, alpha=1, one_tail=True, test=True, norm_fix=True,
                     test_pair=False, y_axis_ran=False, scatter=True, height=False):
    plt.figure(figsize=(p_data.shape[1] / 2, 6))

    c_fr = [np.array([70, 169, 224, 255]) / 255]
    c_cdr = [np.array([239, 92, 85, 255]) / 255]
    c_ot = '#669900'
    c_h100 = '#9966FF'

    if not height:
        height = np.max(np.max(p_data))

    for i in range(p_data.shape[1]):
        if '_fr' in p_data.columns[i] or '_FR' in p_data.columns[i]:
            color = c_fr[0]
            trans = 0.15
        elif '_cdr' in p_data.columns[i] or '_CDR' in p_data.columns[i]:
            color = c_cdr[0]
            trans = - 0.15
        elif 'OT' in p_data.columns[i]:
            color = c_ot
            trans = 0
        elif 'HT' in p_data.columns[i]:
            color = c_h100
            trans = 0

        plot = p_data.iloc[:, i].dropna().astype(float)

        plt.boxplot(plot, widths=0.5, zorder=2, whis=1.5, positions=[i + trans],
                    patch_artist=True, showfliers=False,
                    boxprops={'color': color, 'facecolor': (1, 0, 0, 0), 'linewidth': lw},
                    whiskerprops={'color': color, 'alpha': 0.7, 'linewidth': lw},
                    capprops={'color': color, 'alpha': 0.7, 'linewidth': lw},
                    medianprops={'linestyle': '-', 'color': color, 'alpha': 1, 'linewidth': lw})
        if scatter:
            plt.scatter([i + trans] * len(plot) + (np.random.random(len(plot)) - 0.5) / 3,
                        plot, c=color, s=1, alpha=alpha, edgecolors='none')
        plt.scatter(i + trans, np.mean(plot), marker='+', s=50, color='black', alpha=0.7, zorder=3)

    if test:
        test_info = []
        for i in range(int(p_data.shape[1] / 2)):
            da1 = p_data.iloc[:, i * 2].dropna().astype(float)
            da2 = p_data.iloc[:, i * 2 + 1].dropna().astype(float)
            if test_pair:
                overlap = list(set(list(da1.index)) & set(list(da2.index)))
                overlap.sort()
                da1 = da1.loc[overlap]
                da2 = da2.loc[overlap]
            try:
                xyz = st.mark_significance(plt, da1, da2, i * 2 + 0.15, i * 2 + 0.85, height=height, c='gray',
                                           y_axis_ran=y_axis_ran,
                                           rounding=4, one_tail=one_tail, vertical_ratio=1.5, test_pair=test_pair,
                                           norm_fix=norm_fix)
            except ValueError:
                xyz = [[], [], [], y_axis_ran]

            # 自动获取上一次的yrange，使得全局统一
            test_info += [xyz[2]]
            y_axis_ran = xyz[3]

    ax1 = plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)  # 设置底部坐标轴的粗细
    ax1.spines['left'].set_linewidth(2)  # 设置左边坐标轴的粗细

    # line is a Line2D instance
    for line in ax1.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(2)
    for line in ax1.xaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(2)

    if test:
        return height, test_info
    else:
        return height


def main():
    usage = "usage: python %prog -f file_path -s save_path -i mismatch_path"
    description = "visualization, -f -s -i options are needed"
    parser = argparse.ArgumentParser(prog="%prog 1.0", description=description, usage=usage, add_help=False)
    parser.add_argument("-h", "--help", action="help",
                        help="Show this help message and exit.")

    parser.add_argument("-f", "--file_path", dest="file_path", type=str,
                        help="data file dir-path.")

    parser.add_argument("-s", "--save_path", dest="save_path", type=str,
                        help="save dir-path.")

    parser.add_argument("-p", "--species_path", dest="species_path", type=str,
                        help="species table .txt")

    args = parser.parse_args()
    if not args.file_path:
        parser.print_help()
        sys.exit(1)
    if not args.save_path:
        parser.print_help()
        sys.exit(1)
    if not args.species_path:
        parser.print_help()
        sys.exit(1)

    root_path = args.file_path
    output_path = args.save_path
    species_path = args.species_path

    root_path = st.path_diagnosis(root_path)
    output_path = st.path_diagnosis(output_path)
    st.mkdir(output_path)

    input_path = root_path + 'forward/'
    input_path2 = root_path + 'backward/'

    dpi = 600

    species_order = st.read_file(species_path)
    species_order = species_order.sort_values('order')
    species_order = species_order.sort_values('order2')
    order = st.sort_set_list(list(species_order.iloc[:, 3]))

    test_infomation = []
    for motif in ['C', 'WRC', 'AGC', 'AAC', 'TGC', 'TAC']:
        CDR_plotdata_s = pd.read_csv(input_path + motif + '_CDR_pool.txt', sep='\t')
        CDR_plotdata_as = pd.read_csv(input_path2 + motif + '_CDR_pool.txt', sep='\t')
        FR_plotdata_s = pd.read_csv(input_path + motif + '_FR_pool.txt', sep='\t')
        FR_plotdata_as = pd.read_csv(input_path2 + motif + '_FR_pool.txt', sep='\t')

        columns = CDR_plotdata_s.columns

        CDR_plotdata = pd.concat([CDR_plotdata_s, CDR_plotdata_as])
        CDR_plotdata = CDR_plotdata.reset_index(drop=True)
        CDR_plotdata.columns = [i + '_cdr' for i in columns]
        FR_plotdata = pd.concat([FR_plotdata_s, FR_plotdata_as])
        FR_plotdata = FR_plotdata.reset_index(drop=True)
        FR_plotdata.columns = [i + '_fr' for i in columns]

        plotdata = pd.concat([FR_plotdata, CDR_plotdata], axis=1)

        columns_use = []
        for col in columns:
            columns_use += [col + '_fr', col + '_cdr']
        plotdata = plotdata[columns_use]

        columns_order = []
        for col in order:
            columns_order += [col + '_fr', col + '_cdr']
        plotdata = plotdata.loc[:, columns_order]

        pos_s = np.array(list(range(0, plotdata.shape[1] - 2, 2))) + 0.5  #

        if motif == 'C':
            ymax = 0.15
            stride = 0.05
        elif motif == 'WRC' or motif == 'TAC':
            ymax = 0.3
            stride = 0.1
        elif motif == 'AGC':
            ymax = 0.5
            stride = 0.1
        else:
            ymax = 0.4
            stride = 0.1

        height, test_info = plot_system_pool(plotdata, lw=4, one_tail=False, scatter=False,
                                             test=True, norm_fix=False, height=ymax)
        plt.ylim(-0.02, ymax + 0.03)
        plt.xlim(-1, plotdata.shape[1])
        plt.xticks(pos_s, order, rotation=90, fontsize=15)
        plt.yticks(np.arange(0, ymax + 0.01, stride), fontsize=15)
        plt.title(motif, y=1.1, fontsize=25)

        save_path = output_path + motif + '_mutation_box.pdf'
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)

        height = plot_system_pool(plotdata, lw=4, scatter=False, test=False, height=ymax)
        plt.ylim(-0.02, ymax + 0.03)
        plt.xlim(-1, plotdata.shape[1])
        plt.xticks(pos_s, [''] * len(pos_s), rotation=90, fontsize=15)
        plt.yticks(np.arange(0, ymax + 0.01, stride), [''] * len(np.arange(0, ymax + 0.01, stride)), fontsize=15)

        save_path = output_path + motif + '_mutation_box_unlabel.pdf'
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
        test_infomation += [pd.DataFrame(test_info, index=order)]

    for i, motif in enumerate(['C', 'WRC', 'AGC', 'AAC', 'TGC', 'TAC']):
        test_infomation[i].to_csv(output_path + motif + '_test_info.txt', sep='\t')


if __name__ == "__main__":
    main()
