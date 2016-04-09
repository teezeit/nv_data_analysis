@author: tobias hoelzer
Feb 2016

If you want to try it, clone the repository and run
python run_03.py

Dependencies:
numpy, matplotlib, glob, scipy


QUICKSTART:

Use the python script 'run.py' the following way:
1. import data analysis
import data_analysis
2. Create a sample object initialized with name and path of data folder
sampleNAME = data_analysis.data('NAME','./../DATAFOLDER')
3. Subtract background noise (if necessary)
sampleNAME.set_subtract_noise(subtract_Bool=True/False)
4. Set  i)initial parameters for GAUS FIT
        ii) the run number which is used to fit all x0's
        iii) if sigma should be fitted 
sampleNAME.set_gaus_guess(p_photocurrent = [a,x0,sigma,yo],p_pulsed_laser = [a,x0,sigma,yo],only_a_run_number=NR,fit_sigma=True/False)
5. Run analysis
sampleNAME.analyse_sample()
6. Plot whatever you need
sampleNAME.plot1()
sampleNAME.plot2()
sampleNAME.plot6()


WHAT IS HAPPENING:
1.Sample is initialized with name and directory
2.The directory is parsed an all available data files are saved
3.All the data for single runs is extracted into single variables
4.All the data for single runs is transferred into real physical units
5.For every single run, the peaks are extracted
6.For every single run, the gauss integrals are extracted
  -for this, i need the initial guess parameters, fit the gaussian function y0+a*exp[(x-x0)**2/(2*sigma**2)]
  -you can choose which values to fit
  -the gaussian integral from -inf to inf is then sqrt(PI/a)
7.Then, all the data from the single runs is concatenated into single variables for the sample



Additionally you can compare samples using the Experiment class:

1. Initialize an experiment object with data objects SAMPLE1 and SAMPLE2
Experiment = data_analysis.experiment(SAMPLE1, SAMPLE2)
2. Plot both Samples
Experiment.plot_experiment()







