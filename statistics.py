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



class cox_engine(object):

    def __init__(self, data, time, censoring):

        self.order_of_time = np.argsort(time,axis=0).flatten()
        self.num_times, self.num_vars = data.shape

        self.data = data[self.order_of_time]
        self.times = time[self.order_of_time]
        self.censoring = censoring[self.order_of_time].astype(int)

        self.b = np.zeros((1,len(self.data[0,:])))

    def get_exponential_sum(self):
        self.exp_sum = np.exp(np.sum(self.b * self.data,axis=1,keepdims=1))

    def get_g(self):
        self.g = np.sum(self.exp_sum,axis=0) - np.cumsum(self.exp_sum,axis=0) + self.exp_sum

    def get_h(self):
        inner_sum = self.data * self.exp_sum
        self.h = np.sum(inner_sum,axis=0) - np.cumsum(inner_sum,axis=0) + inner_sum

    def get_u(self):
        self.u = np.sum(self.censoring * (self.data - self.h/self.g),axis=0,keepdims=1)

    def get_a(self):
        pass
        inner_sum = np.transpose(self.data * np.ones((self.num_vars,self.num_times,self.num_vars)) * self.exp_sum)
        inner_sum = self.data * np.ones((1,self.num_times,self.num_vars)) * inner_sum
        misordered_sum = np.cumsum(inner_sum[:,::-1,:],axis=1)[:,::-1,:]

        self.a = np.ones((self.num_vars,self.num_times,self.num_vars))
        for _ in range(self.num_times):
            self.a[::1,_,:] = misordered_sum[::1,_,:]

    def get_i(self):
        inner_sum = self.data * np.transpose(self.data * np.ones((self.num_vars,self.num_times,self.num_vars)) / self.g)
        hg = np.ones((self.num_vars,self.num_times,self.num_vars))
        for _ in range(self.num_times):
            hg[::1,_,:] = inner_sum[::1,_,:]
        inner_sum = self.a - hg
        inner_sum *= (np.ones((self.num_times,self.num_vars)) * self.censoring/self.g)
        self.i_inv = np.linalg.inv(-1*np.sum(inner_sum,1))

    def calc_LL(self):
        return np.sum(self.censoring * (np.sum(self.b * self.data,1) - np.log(self.g)))

    def get_all(self):
        self.get_exponential_sum()
        self.get_g()
        self.get_h()
        self.get_u()
        self.get_a()
        self.get_i()

        self.b = (self.b.transpose() - np.dot(self.i_inv, self.u.transpose())).transpose()

    def iterate(self):
        self.get_all()
        last_step = self.calc_LL()
        current_step = last_step + 1
        while np.abs(last_step - current_step) > 1e-8:
            last_step = 1.*current_step
            print("Iterating!")
            self.get_all()
            current_step = self.calc_LL()

