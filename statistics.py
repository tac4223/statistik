# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 17:31:13 2015

@author: xerol
"""
import numpy as np

class analysis(object):

    def __init__(self, filename, headersize):

        self.raw_data = np.loadtxt(filename,skiprows = headersize)
        covars = "Alter	Geschlecht	T N GTVVol FDGVol	TBR_0	V16_0	\
        TBR_1	V16_1	TBR_2	V16_2	TBR_5	V16_5".split()

        self.subject_count = len(self.raw_data[:,0])

        self.covar_names = {}
        for num in range(14):
            self.covar_names[num + 3] = covars[num]


    def input_chain(self,prompt,expected_type):
        query = raw_input(prompt)
        try:
            expected_type(query)
            return expected_type(query)
        except:
            print("{0} ist keine gültige Eingabe, bitte wiederholen".\
                format(query))
            self.input_chain(prompt, expected_type)

    def pick_cohorts(self):
        print("Bitte gewünschte Kohorte bzw. Kohortenkombination wählen!")
        cohort = self.input_chain("0: Experimentalkohorte\n1: Validierung\n2: "\
            "kombiniert\nAuswahl: ", int)
        if cohort in [0,1]:
            self.cohort_mask = (self.raw_data[:,0] == cohort+1)
        elif cohort == 2:
            self.cohort_mask = np.ones(self.subject_count,dtype=bool)
        else:
            print("Uneindeutige Eingabe, bitte wiederholen!\n")
            self.pick_cohorts()
        self.test_times = np.reshape(self.raw_data[:,1][self.cohort_mask],
             (-1,1))
        self.test_censored = np.reshape(self.raw_data[:,2][self.cohort_mask],
            (-1,1))

    def pick_vars(self):
        go_on = 1
        variables = []
        while go_on:
            print("Bitte gewünschte Variable auswählen:\n")
            for num in range(14):
                print("{0}: {1}".format(self.covar_names[num+3],num))

            usr_input = self.input_chain("Bitte Zahl (0-13) eingeben: ", int)
            if usr_input <= 13:
                variables.append(usr_input + 3)
            else:
                print("Nur Zahlen bis 13 eingeben!")

            keep_going = self.input_chain("Weitere Variable auswählen? (j/n)",
                                          str)

            if keep_going in ["j","n"]:
                go_on = (keep_going == "j")
            else:
                print("Unklare Eingabe, springe zurück zur Variablenauswahl...")

        print("\nEventuell doppelt eingegebene Variablen werden nur einmal "\
        "berücksichtigt.")
        self.test_variables = [self.covar_names[num] for num in
            np.unique(variables)]
        self.test_data = self.raw_data[:,np.unique(variables)][self.cohort_mask]




datacore = analysis("Beleg_Biostatistik_2015_Daten.txt",7)
datacore.pick_cohorts()
datacore.pick_vars()
for _ in [datacore.test_censored,datacore.test_times,datacore.test_variables, datacore.test_data]:
    print(_)