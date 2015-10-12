# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 17:31:13 2015

@author: xerol
"""
import numpy as np
import scipy.stats as sp

class analysis(object):
    """
    Klasse um Daten für statistische Analysen vorzubereiten. Liest CSV-Dateien
    ein, stellt sicher dass Nutzereingaben sinnvoll sind und ruft schließlich
    die eigenliche statistische Analyse auf.

    Funktionen:
    fit_cox: Ruft die Klasse mit Cox-Regression auf.
    fit_loop: Gibt dem Nutzer nach jedem Analysedurchlauf die Möglichkeit,
        eine weitere Analyse anzustoßen oder das Programm zu beenden.
    fit_statistics: Berechnet die in der Aufgabenstellung geforderten Parameter
        auf Basis der Fitparameter b_j und gibt sie (mehr oder weniger) über-
        sichtlich aus.
    input_chain: Hauptfunktion für Userinput. Stellt sicher dass gewünschte
        Eingabe getätigt wird.
    pick_cohorts: Funktion die den Nutzer auswählen lässt welche Kohorten
        betrachtet werden sollen.
    pick_vars: Auswahl der zu fittenden Kovariablen.
    stat_to_file: Funktion um Daten aus fit_statistics nicht in Konsole auszu-
        geben sondern direkt in Logfile zu schreiben.

    Instanzvariablen:
    cohort_mask: Maske die Daten basierend auf der gewählten Kohorte auswählt.
    covar_names: Dictionary das einen Zusammenhang zwischen Spalten in
        Datenarray und den zugehörigen Namen herstellt.
    hazard_confidence: Konfidenzintervalle der berechneten Hazard Ratios.
    Hazard Ratio: Genau das.
    p_val: p-Werte der einzelnen Fitparameter, ermöglicht Aussage über
        Signifikanz.
    raw_data: Array mit allen Daten, hieraus wird Datenarray für konkrete
        Analyse aufgefüllt.
    significance: Nutzereingabe, Signifikanzniveau für Chi²-Test.
    stabw: Standardabweichung der Fitparameter.
    subject_count: Anzahl der Einträge je Kovariable.
    test_censored: Enthält Zensierungen für die Regression
    test_data: Array mit Daten der einzelnen Kovariablen für die Regression.
    test_times: Array mit Zeiten für die Regression
    test_variables: Namen der getesteten Kovariablen. Sonst hat man ja
        hinterher keine Ahnung mehr.
    wald: Wald-Statistik der Fitparameter.

    Instanzen:
    fit: cox_engine-Instanz
    """

    def __init__(self, filename, headersize):
        """
        Initialisierung. Filename ist CSV-Tabelle, headersize die Anzahl der
        Nicht-numerischen Zeilen vor den eigentlichen Daten.
        """
        self.raw_data = np.loadtxt(filename,skiprows = headersize)
        covars = "Alter	Geschlecht	T	N	GTVVol	FDGVol	TBR_0\
        V16_0	TBR_1	V16_1	TBR_2	V16_2	TBR_5	V16_5\
        cutoff_tbr2".split()
        self.num_of_covars = len(covars)

        self.subject_count = len(self.raw_data[:,0])

        self.covar_names = {}
        for num in range(len(covars)):
            self.covar_names[num + 3] = covars[num]


    def fit_loop(self):
        """
        Ausgangspunkt, ruft alle relevanten Funktionen auf und gibt dann Mög-
        lichkeit zum erneuten Ausführen.
        """
        go_on = 1
        while go_on:
            self.pick_cohorts()
            self.pick_vars()
            self.fit_cox()
            self.fit_statistics()
            self.print_statistic()
            self.stat_to_file()
            go_on = self.input_chain("Neue Analyse starten? (j/n): ",str) == "j"

    def mass_analysis(self, number):
        """
        number: Index der gewünschten klinischen Variable zur Erstellung aller
            Hypoxie-Modelle, entsprechend self.pick_vars()

        Reine Backend-Funktion um die massenhafte Erstellung von Modellen mit
        je einer Hypoxievariable zu vereinfachen. Fire & Forget, die Variable
        die per number übergeben wird, wird einmal komplett durchgehauen.
        """
        self.cohort_mask = np.ones(self.subject_count,dtype=bool)
        self.test_times = np.reshape(self.raw_data[:,1][self.cohort_mask],
             (-1,1))
        self.test_censored = np.reshape(self.raw_data[:,2][self.cohort_mask],
            (-1,1))
        for _ in range(6,self.num_of_covars-2):
            variables = [number+3, _+3]
            self.test_variables = [self.covar_names[num] for num in variables]
            self.test_data = self.raw_data[:,np.unique(variables)]\
                [self.cohort_mask]
            self.fit_cox()
            self.fit_statistics(0.05)
            self.stat_to_file()


    def input_chain(self,prompt,expected_type):
        """
        prompt: String der als Prompt für den Nutzer angezeigt werden soll.
        expected_type: Welcher Variablentyp als Eingabe akzeptiert wird.

        Stellt standardisierte Möglichkeit für Nutzereingaben zur Verfügung.
        """
        query = raw_input(prompt)
        try:
            return expected_type(query)
        except:
            print("{0} ist keine gültige Eingabe, bitte wiederholen".\
                format(query))
            self.input_chain(prompt, expected_type)


    def pick_cohorts(self):
        """
        Ermöglicht die Nutzerauswahl der zu betrachtenden Kohorten.
        """
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
        """
        Nutzerauswahl der zu berücksichtigenden Kovariablen.
        """
        go_on = 1
        variables = []
        while go_on:
            print("Bitte gewünschte Variable auswählen:\n")
            for num in range(self.num_of_covars):
                print("{0}: {1}".format(self.covar_names[num+3],num))

            usr_input = self.input_chain("Bitte Zahl (0-{0}) eingeben: ".
                format(self.num_of_covars-1), int)
            if usr_input <= self.num_of_covars-1:
                variables.append(usr_input + 3)
            else:
                print("Nur Zahlen bis {0} eingeben!".format(
                    self.num_of_covars-1))

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


    def fit_cox(self):
        """
        Tut im Prinzip nichts außer die Cox-Engine zu initialisieren und
        anschließend die Iteration anzustoßen.
        """
        self.fit = cox_engine(self.test_data, self.test_times,
              self.test_censored)
        self.fit.iterate()


    def fit_statistics(self, significance=None):
        """
        Berechnet alle in der Aufgabenstellung geforderten Werte für die
        gefundenen Fitparameter. Sollte eigentlich nichts unverständliches
        dabei sein.
        """
        if significance == None:
            self.significance = self.input_chain("Bitte gewünschtes "\
                "Signifikanzniveau in Prozent eingeben: ", float)/100.
        else: self.significance = significance
        self.stabw = np.sqrt(np.diagonal(-1*self.fit.i_inv))
        self.wald = (self.fit.b/self.stabw)**2
        self.p_val = 1 - sp.chi2.cdf(self.wald,1)
        self.hazard_ratio = np.exp(self.fit.b)
        self.hazard_confidence = self.hazard_ratio - np.exp(sp.norm.ppf(1 -
            self.significance/2) * self.stabw)
        self.hazard_confidence = np.array([(self.hazard_ratio -
            self.hazard_confidence).flatten(),
            (self.hazard_ratio + self.hazard_confidence).flatten()])


    def print_statistic(self):
        """
        Gibt die in fit_statistics berechneten Daten in Konsole aus.
        """
        print("Gefittete Parameter:\n\n")
        for param in range(self.fit.num_vars):
            print("\n{0}\n-----\nb: {1}\nSt.abw.: {2}\nHazard-Ratio: {3}\n"\
            "Hazard-KI: {4}\nWald-Statistik: {5}\np-Wert: {6}".format(
            self.test_variables[param],self.fit.b[0,param],self.stabw[param],
            self.hazard_ratio[0,param],self.hazard_confidence[:,param],
            self.wald[0,param],self.p_val[0,param]))
            if self.p_val[0,param] < 0.05:
                print("Signifikanter Einfluss!")
            else:
                print("Kein signifikanter Einfluss!")


    def stat_to_file(self, filename="results.txt"):
        """
        filename: string, Name der Datei die angelegt werden soll. Falls schon
            vorhanden, werden neue Ausgaben angehängt.

        Schreibt die Daten in Logfile statt sie auf Konsole auszugeben.
        """
        with open(filename,"a") as logfile:
            logfile.write("\n\n################\n{}\n################\n".\
                format(self.test_variables))
            for param in range(self.fit.num_vars):
                logfile.write("\n-----\n{0}\n-----\nb: {1}\nSt.abw.: {2}\n"\
                "Hazard-Ratio: {3}\n"\
                "Hazard-KI: {4}\nWald-Statistik: {5}\np-Wert: {6}".format(
                self.test_variables[param],self.fit.b[0,param],
                self.stabw[param],
                self.hazard_ratio[0,param],self.hazard_confidence[:,param],
                self.wald[0,param],self.p_val[0,param]))


class cox_engine(object):
    """
    Named after Zefram Coxrane, who built humanity's first regression-capable
    vessel out of an old Titan II nuclear p-value.

    Führt (hoffentlich) eine Cox-Regression aus, entsprechend der Anleitung in
    der Aufgabenstellung.

    Funktionen:
    calc_LL: Berechnet die Log-Likelihood für derzeitige b_j
    get_a: Berechnet nxnxn-Matrix A, die für den Fit benötigt wird.
    get_all: Ruft alle "get" Funktionen in korrekter Reihenfolge auf,
        aktualisiert b nach jedem Durchlauf.
    get_exponential_sum: Die hier berechnete Summe wird immer wieder benötigt,
        daher erschien es sinnvoll die Berechnung bequem auszulagern.
    get_g: Berechnet g entsprechend Anleitung.
    get_h: Berechnet h.
    get_i: Berechnet die Matrix I, die für die iterative Annäherung von b
        nötig ist.
    get_u: Berechnet den Vektor u, der ebenfalls direkt in die Näherung für b
        eingeht.
    iterate: Solange nähern bis Schwelle unterschritten.

    Variablen:
    a: Matrix A, kubisch, siehe Anleitung (8)
    b: Die aktuellen Werte für b, wird mit 0 initialisiert.
    censoring: Zensierungen für alle Datenpunkte.
    data: Zu fittende Daten.
    exp_sum: Exponential von Summe über Achse 1 von b*data
    g: g-Vektor, siehe Anleitung (6).
    h: h, siehe Anleitung (7)
    i_inv: Wie Anleitung (10), aber schon invertiert.
    num_times: Anzahl der Datenpunkte pro Kovariable
    num_vars: Anzahl der Kovariablen.
    order_of_time: Sortierindizes um Daten in Reihenfolge aufsteigender Zeiten
        zu bringen.
    times: Zeitpunkte für jede Kovariable.
    u: Vektor u, siehe Anleitung (9)
    """

    def __init__(self, data, time, censoring):
        """
        data: np-Array mit den Daten innerhalb einer Kovariable.
        time: Array mit zugehörigen Zeiten.
        censoring: Array mit Zensierungen.

        Leistet eigentlich nur die Initialisierung der Objektinstanz, ansonsten
        geschieht nicht viel spannendes.
        """
        self.order_of_time = np.argsort(time,axis=0).flatten()
        self.num_times, self.num_vars = data.shape

        self.data = data[self.order_of_time]
        self.times = time[self.order_of_time]
        self.censoring = censoring[self.order_of_time].astype(int)

        self.b = np.zeros((1,len(self.data[0,:])))

    def get_exponential_sum(self):
        """
        Summe wird immer wieder benötigt, lässt sich gut auslagern.
        """
        self.exp_sum = np.exp(np.sum(self.b * self.data,axis=1,keepdims=1))


    def get_g(self):
        """
        g berechnen, siehe Anleitung.
        """
        self.g = np.sum(self.exp_sum,axis=0) - np.cumsum(self.exp_sum,axis=0)\
            + self.exp_sum


    def get_h(self):
        """
        h berechnen, entsprechend Anleitung.
        """
        inner_sum = self.data * self.exp_sum
        self.h = np.sum(inner_sum,axis=0) - np.cumsum(inner_sum,axis=0) + \
            inner_sum


    def get_u(self):
        """
        u berechnen, auch eher straight forward.
        """
        self.u = np.sum(self.censoring * (self.data - self.h/self.g),axis=0,
            keepdims=1)


    def get_a(self):
        """
        Berechnet A... glaube ich. Dreidimensionale Matrixwürfel verwirren mich,
        Numpy-Broadcasting verwirrt mich und das Ergebnis könnte völliger
        Unfug sein. Use with care.

        inner_sum: x_mj*x_mk*exp... entsprechend Anleitung
        misordered_sum: Ergebnis von völlig verstrahltem Broadcasting.
            Vermutlich käme man schöner an das richtige Ergebnis, sofern
            mein Ergebnis überhaupt korrekt ist. Mir fiel nur nichts mehr ein.
        """
        inner_sum = np.transpose(self.data * np.ones((self.num_vars,
              self.num_times,self.num_vars)) * self.exp_sum)
        inner_sum = self.data * np.ones((1,self.num_times,self.num_vars)) *\
            inner_sum
        misordered_sum = np.cumsum(inner_sum[:,::-1,:],axis=1)[:,::-1,:]

        self.a = np.ones((self.num_vars,self.num_times,self.num_vars))
        for _ in range(self.num_times):
            self.a[:,_,:] = misordered_sum[:,_,:]


    def get_i(self):
        """
        Falls A korrekt ist, sollte hier nicht viel schief gehen. Aber wieder
        hatte ich mehrere große Knoten im Hirn beim Versuch, das vernünftig
        hinzukriegen.

        inner_sum ist der Broadcast von h_ij*h_ik auf den entsprechenden
            Matrixwürfel, dabei kommen allerdings die Elemente in der falschen
            Reihenfolge zutage (alle ersten Zeilen, alle zweiten Zeilen...)
        hg der identische Würfel wie inner_sum, aber korrekt sortiert. Geht
            sicher schöner.
        """
        inner_sum = self.data * np.transpose(self.data * np.ones(
            (self.num_vars,self.num_times,self.num_vars)) / self.g)
        hg = np.ones((self.num_vars,self.num_times,self.num_vars))
        for _ in range(self.num_times):
            hg[:,_,:] = inner_sum[:,_,:]
        inner_sum = self.a - hg
        inner_sum *= (np.ones((self.num_times,self.num_vars)) * \
            self.censoring/self.g)
        self.i_inv = np.linalg.inv(-1*np.sum(inner_sum,1))


    def calc_LL(self):
        """
        Gibt die LL-Funktion für jeden Wert in self.b zurück.
        """
        return np.sum(self.censoring * (np.sum(self.b * self.data,1) -
            np.log(self.g)))


    def get_all(self):
        """
        Ruft alles auf was man zum Fitten so braucht.
        Berechnet alle Parameter, aktualisiert zu guter letzt die b_j
        """
        self.get_exponential_sum()
        self.get_g()
        self.get_h()
        self.get_u()
        self.get_a()
        self.get_i()

        self.b = (self.b.transpose() - np.dot(self.i_inv, self.u.transpose()))\
            .transpose()


    def iterate(self, threshold=1e-8):
        """
        Threshold: Mindeständerung in der LL pro Iterationsschritt um den Fit
        noch nicht abzubrechen.

        last_step: Wert der LL Funktion aus dem letzten Schritt.
        current_step: Entsprechender derzeitiger Wert.
        """
        self.get_all()
        last_step = self.calc_LL()
        current_step = last_step + 1
        while np.abs(last_step - current_step) > threshold:
            last_step = 1.*current_step
            self.get_all()
            current_step = self.calc_LL()

datacron = analysis("Beleg_Biostatistik_2015_Daten_tbrcutoff.txt",8)
datacron.fit_loop()