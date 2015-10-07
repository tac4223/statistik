# -*- coding: utf-8 -*-
"""
Created on Tue Oct 06 21:21:12 2015

@author: xerol
"""
import statistics

datacore = analysis("Beleg_Biostatistik_2015_Daten.txt",7)
datacore.pick_cohorts()
datacore.pick_vars()
for _ in [datacore.test_censored,datacore.test_times,datacore.test_variables, datacore.test_data]:
    print(_)