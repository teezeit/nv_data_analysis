#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Feb 2016
This file uses the data_analysis.py framework to Analyse the Experimental Data
of two Diamond Samples, one with and one withou NV-Centers.


@author: tobias hoelzer
"""

#imports
import data_analysis




#Create data_analysis object [we initialize with the sample name and the directory where the files are stored]
sample1 = data_analysis.data('Sample1: Jeson CVD NV','./../data/03-09-16/jeson_nv/')
sample4 = data_analysis.data('Sample4:  HPHT Non nv','./../data/03-04-16/2dnonnv/')


#subtract noise from data (default is None)
sample1.set_subtract_noise(subtract_Bool=True)
sample4.set_subtract_noise(subtract_Bool=True)

#Set the initial guess parameters for the gaus fitting functions:
#photodiode,
#gaus(x,a,x0,sigma,y0)-> P = [a,x0,sigma,y0], only_a_run_number determins which run is used to fit the other 3 parameters
sample1.set_gaus_guess(p_photocurrent=[0.00005,500., 50.,-0.000024],p_pulsedpower =[100.,500., 50.,0.],only_a_run_number=3,fit_sigma = False)
sample4.set_gaus_guess(p_photocurrent=[0.00005,500., 50.,-0.000024],p_pulsedpower =[100.,500., 50.,0.],only_a_run_number=3,fit_sigma = False)

sample1.analyse_sample()
sample4.analyse_sample()

#sample4.gaus_fits()
#sample4.plot()

#sample4.plot2()
#sample1.plot3()
#sample1.plot4()
#sample1.plot5()
#sample4.plot6()



Experiment =data_analysis.experiment(sample4,sample1)
#print Experiment.samples

Experiment.plot_experiment()
