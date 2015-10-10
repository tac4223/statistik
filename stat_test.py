# -*- coding: utf-8 -*-
"""
Created on Tue Oct 06 21:21:12 2015

@author: xerol
"""
import statistics as st

datacore = st.analysis("Beleg_Biostatistik_2015_Daten.txt",7)
datacore.pick_cohorts()
datacore.pick_vars()

cox = st.cox_engine(datacore.test_data, datacore.test_times, datacore.test_censored)
