# -*- coding: UTF-8 -*-
# @Author: Zhang Senxin

"""
unclassified code
now only comtain NPS summary
"""

import pandas as pd


def summary_system(pathout, summary):
    isExists = os.path.exists(pathout + 'summary.txt')
    if isExists:
        summary0 = []
        for line in open(pathout + 'summary.txt'):
            summary0 += line.strip('\n').split('\t')
        if not (summary[0] in summary0):
            summary0 += summary
    else:
        summary0 = summary
    summary0 = pd.DataFrame(summary0)
    summary0.to_csv(pathout + "summary.txt", sep='\t', index=False, header=False)
