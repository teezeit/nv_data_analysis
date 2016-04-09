#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Feb 2016
This is the  data_analysis.py framework

1. The master class is Experiment wich contains
2. sample objects wich contains
3. data objects wich contains a list of 
4. single run objects


The methods from the data class mostly loop throught the list of runs and call the respective functions there.
The rest is readin and plotting
every textfile in the sample folder is another single run




How to use data_analysis:
1. Create Sample Object, initialize with the Sample name and the folder path were all the experimental runs are
    Sample01=data('Sample One','./')
2. Set the initial parameters for the Gaus fit used to calculate the integral of the photocurrent and the pulsedpower,
    define which run number is used to fit the initial gaus curve that is used, define if only gaus height or also sigma is fitted
        Sample01.set_gaus_guess(p_photocurrent=[0.000005,600., 50.,-0.000024],p_pulsedpower =[100.,600., 50.,0.], only_a_run_nuber=3,fit_sigma=False)
3. Run Sample.analyse_sample()
    Sample01.analyse_sample()
4. Plot whichever data you want to see per sample
    Sample01.plot4()
    Sample02.plot5()

5. You can initialize an Experiment object with multiple sample objects and plot them together Experiment.plot_experiment
    Experiment =experiment(sample01)
    Experiment.plot_experiment()


@author: tobias hoelzer
"""

import numpy as np                             #to play with data
import matplotlib.pyplot as pt                 #to plot data
import matplotlib as mpl                     #used to change colorsequence of plots
from mpl_toolkits.mplot3d import Axes3D     #plot 3D plots
import glob                                    #to search for textfiles in folder
from scipy.optimize import curve_fit


#optional imports:
#seaborn is for nice looking plots (but if activated, you cannot change the fontsize, etc.)
#import seaborn as sns
#matplotlib.cm is if you want to import colormaps
#import matplotlib.cm as cm



#This is an Experiment class. It contains several sample classes.
#Each sample is a different diamond (in our case 1 NV and 1 non NV)
class experiment(object):
    
    def __init__(self,*args):
        self.samples = args

    #This method 3d plots the cw power VS. Pulsed laser power peak VS. Photocurrent Peak
    #for all samples contained in the 'experiment object'.
    def plot_experiment(self):

        #create figure
        fig = pt.figure()
        #create subplot
        ax = fig.add_subplot(111,projection = '3d')
        #create legend list
        ll = []
        
        #loop through samples to input data
        for (j,sample) in enumerate(self.samples):
            x = sample.S_cw_power
            y = sample.S_pulsedpower_integral
            z = sample.S_photocurrent_integral
            ll.append(sample.name)

            if j==0:
                ax.scatter(x,y,z,c='b',marker='o')
            elif j==1:
                ax.scatter(x,y,z,c='r',marker ='^')
            elif j==2:
                ax.scatter(x,y,z,c='g',marker ='x')
            ax.grid()

        #Set Labels
        ax.set_xlabel('cw power [\mu W]')
        ax.set_ylabel('Laser Power per Pulse [mJ]')
        ax.set_zlabel('Photocurrent Integral [nC]')
        ax.set_title('Sample Comparison')
        ax.legend(ll)
        #pt.rc('font',size=30)


        
        ##2D plot###################################################
        

        fig2 = pt.figure()
        #loop through samples to input data
        for (j,sample) in enumerate(self.samples):
            
            subplot_i = str(len(self.samples)) +'1' + str(j+1)
            #add second subplot
            ax2 = fig2.add_subplot(int(subplot_i))
            #create legend list
            ll2 = []
            if     j == 0:
                cm = pt.get_cmap('Blues')
                m = 'o'
            elif j == 1:
                cm = pt.get_cmap('Oranges')
                m = '^'
            elif j == 2:
                cm = pt.get_cmap('copper')
                m = 'x'
                
            #Set colorcycle
            ax2.set_color_cycle([cm(1.*i/len(sample.runs)) for i in range(len(sample.runs))])
            
            #iterate through all runs
            for run in sample.runs:
                #get data
                #we want to plot them in order of their pulsedpower so sort (imortant along axis0)!
                indexorder = run.pulsedpower_integral.argsort(axis=0)
                #for some reason the dimensions get fucked, squeeze removes unnecessary dimensions.
                x2 = np.squeeze(run.pulsedpower_integral[indexorder])
                y2 = np.squeeze(run.photocurrent_integral[indexorder])
                #save legend entry
                ll2.append('CW [uW] = ' + str(run.cw_power) )
                #plot
                ax2.plot(x2,y2,linestyle ='-',linewidth = 3, marker = m)
                
        
            ax2.grid()
            #Set Labels
            ax2.set_xlabel('Laser Power per Pulse [mJ]')
            ax2.set_ylabel('Photocurrent Integral [nC]')
            ax2.set_title('Sample: ' +  sample.name)
            ax2.legend(ll2,prop={'size':8})
        
        pt.show()
        




#This is a sample class. a Sample is an NV or an nonNV data set
class sample(experiment):

    #Initialize
    def __init__(self,samplename,directory):
        super(sample,self).__init__()
        #Name of the sample
        self.name = samplename

        #the directory of where the data is located
        self.dir = directory

        

#this is a data class. a data class contains
#1) a list of data files and 2) a list of single_run objects 
class data(sample):
    
    #Initialize
    def __init__(self,*args):
        #call parent class and initialize data with sample name and directory
        super(data,self).__init__(*args)

        #list of filenames in this directory. each file contains a different CW run
        self.filenames = glob.glob(self.dir + '*.txt')
        self.filenames.sort()

        #list of run objects. each run has a different CW power
        #and contains raw time, photocurrent data and photodiode data
        self.runs = None

        #initial guess for gaus fit a*exp[(x-x0)^2/(2sigma^2)]+y0
        self.gaus_guess =[[0.000005,600., 50.,-0.000024],[100.,600., 50.,0.]]
        #initial guess for fit a*exp[(x-x0)^2/(2sigma^2)]+y0 has been calculated or guessed?
        #if True, then we take x0,sigma and y0 for granted
        self.fit_only_a = None
        #this dude determines which run is used to calculate the parameters x0,sigma,y0 for gaus fit
        self.only_a_run_number = None
        #this determines if sigma is fitted to gaus or only a
        self.fit_sigma = None
        
        #remove noise from data [only if you measured noise]        
        self.subtract_noise_Bool = None;
        
        #concatenation of measurement data [0] will later be removed. we need it to use concatenate
        self.S_photocurrent_max         = [0]
        self.S_pulsedpower_max        = [0]
        self.S_photocurrent_max_avg     = [0]
        self.S_pulsedpower_max_avg        = [0]
        self.S_cw_power                 = [0]
        self.S_photocurrent_integral     = [0]
        self.S_pulsedpower_integral     = [0]

    # THIS BADBOY  METHOD IS CALLING ALL THE HEAVY LIFTING
    def analyse_sample(self):
        #Create the list of raw data objects
        #read_all_data_ fills object list 'runs' with the 'single_run' objects which contain the raw data
        self.read_all_data()
        #if we measured noise before we want to subtract it from the data
        if self.subtract_noise_Bool == True: self.subtract_noise_from_runs()
        #now extract time, photocurrent, photodiode  peaks from the raw data
        self.calc_real_units()
        if self.fit_only_a == True: self.calc_gaus_guess_for_photocurrent(self.only_a_run_number)
        self.calc_peaks()
        self.calc_integrals(self.gaus_guess,self.fit_only_a,self.fit_sigma)
        # and concatenate all runs into single variables
        self.concatenate_runs()

    # reads in all datatextfiles found in the sample folder, reads them in and extracts the data
    # also reads in the CW power which is only one file for all runs
    def read_all_data(self):

        #read in CW data [mW]
        cw_power_load = np.loadtxt(self.dir + 'cw')
        #cut of first column containing the run number
        cw_power_load = cw_power_load[:,1]


        #create list of run objects 
        self.runs = [single_run(i) for i in self.filenames]

        #get the data from the raw data for all CW runs
        
        for (j,run) in enumerate(self.runs):
            run.extract_data()
            run.set_cw(cw_power_load[j])
    
    
    #Removes noise
    def subtract_noise_from_runs(self):
        
        #read in noise
        noise_raw = np.loadtxt(self.dir + 'noise')
        #cut out the pc data: get every second row starting with the second. remember numbering [0,1,2,3]
        photocurrent_noise = noise_raw[1::2]
        #get mean noise
        photocurrent_noise_mean = np.mean(photocurrent_noise,0)        
                
        
        #remove noise from every
        for run in self.runs:
            run.subtract_noise(photocurrent_noise_mean)
        
    
    #Sets initial gaus parameters
    def set_gaus_guess(self,p_photocurrent,p_pulsedpower,only_a_run_number,fit_sigma):
        self.gaus_guess = [p_photocurrent,p_pulsedpower]
    
        #remember that inital fit paremeters shall be calculated so we only fit the height of the gaus
        if only_a_run_number is not None:
            self.only_a_run_number = only_a_run_number
            self.fit_only_a = True
            
        
        if fit_sigma == True:
            self.fit_sigma = True
        
    #calculates the x0,sigma and y0 position of the GaussFit
    def calc_gaus_guess_for_photocurrent(self,run_number):
        #get the run
        run = self.runs[run_number]
        #define x and y values
        x = run.time
        y = np.mean(run.photocurrent_real,0)
        
        #get initial guess for photocurrent
        p_photocurrent =self.gaus_guess[0]
        
        #fit gaus curve
        popt,pcov = curve_fit(gaus,x,y,p0=p_photocurrent)
        
        #update initial gausfit parameters
        self.gaus_guess[0]    = popt
        
        #plot fit
        if False:
            print 'gauss guess calc:', popt
            pt.close()
            fig = pt.figure()
            ax = fig.add_subplot(111)
            ym= gaus(x,popt[0],popt[1],popt[2],popt[3])
            ax.plot(x,ym,c='r',label = 'best fit')
            ax.scatter(x,y,c='b',label='data photocurrent')
            ax.legend()
            fig.show()
            raw_input('gaus guess calc. press enter to continue...')
            pt.close()
        
        
        
        
                    
        
    #Set subtract noise parameter
    def set_subtract_noise(self,subtract_Bool):
        self.subtract_noise_Bool = subtract_Bool
            


    #plots photocurrent or photodiode of run defined in 'number'
    #yvalue = {'photocurrent', 'photodiode', 'photocurrent_real','photodiode_real'}
    def plot_xy(self,y_value,run_number):
        self.runs[run_number].plot_xy(y_value)

    #for all runs, transforms input data into real SI measurement units 
    def calc_real_units(self):
        for run in self.runs:
            run.return_real_units()

    #for all runs, get the peaks 
    def calc_peaks(self):
        for run in self.runs:
            run.return_peaks()
            
    #for all runs, get the integrals
    def calc_integrals(self,gaus_guess,fit_only_a,fit_sigma):
        for run in self.runs:
            run.return_integrals(gaus_guess,fit_only_a,fit_sigma)

    #concatenate all runs into single data lists
    def concatenate_runs(self):
        
        #go trough all runs and concatenates the data into single arrays
        for run in self.runs:
            
            #those are the peaks of the runs
            self.S_photocurrent_max         = np.append(self.S_photocurrent_max,run.photocurrent_max)
            self.S_pulsedpower_max         = np.append(self.S_pulsedpower_max,run.pulsedpower_max)
            
            #those are the integrals of the runs
            self.S_photocurrent_integral    = np.append(self.S_photocurrent_integral,run.photocurrent_integral)
            self.S_pulsedpower_integral        = np.append(self.S_pulsedpower_integral,run.pulsedpower_integral)

            #this is the  average per sun
            self.S_photocurrent_max_avg     = np.append(self.S_photocurrent_max_avg,np.repeat(np.average(run.photocurrent_max),len(run.photocurrent_max)))
            self.S_pulsedpower_max_avg         = np.append(self.S_pulsedpower_max_avg,np.repeat(np.average(run.pulsedpower_max),len(run.pulsedpower_max)))
            
            #this is the CW power
            self.S_cw_power = np.append(self.S_cw_power,np.repeat(run.cw_power,len(run.pulsedpower_max)))
            #self.photocurrent_max = [self.photocurrent_max,run.photocurrent_max]
            #self.pulsedpower_max  = [self.pulsedpower_max,run.pulsedpower_max]
            
            #don't forget to delete the [0] in the first entry.
            self.S_photocurrent_max         = np.delete(self.S_photocurrent_max,0,0)
            self.S_pulsedpower_max         = np.delete(self.S_pulsedpower_max,0,0)
            self.S_photocurrent_integral    = np.delete(self.S_photocurrent_integral,0,0)
            self.S_pulsedpower_integral        = np.delete(self.S_pulsedpower_integral,0,0)
            self.S_photocurrent_max_avg     = np.delete(self.S_photocurrent_max_avg,0,0)
            self.S_pulsedpower_max_avg        = np.delete(self.S_pulsedpower_max_avg,0,0)
            self.S_cw_power                 = np.delete(self.S_cw_power,0,0)





    #Plot: max photocurrent & pulsedpower over time
    def plot(self):
        fig, ax1 = pt.subplots()
        ax1.plot(self.S_photocurrent_integral,'b',self.S_photocurrent_max_avg,'bs')
        ax2 = ax1.twinx()
        ax2.plot(self.S_pulsedpower_integral,'r',self.S_pulsedpower_max_avg,'rs')
        ax1.set_xlabel('run')
        ax1.set_ylabel('Photocurrent')
        ax2.set_ylabel('Pulsed Laser Power')
        ax1.legend(['Photocurrent Peak', 'Photocurrent Mean'],'upper center')
        ax2.legend(['Pulsed Laser Power Max','Pulsed Laser Power Mean'],'upper right')
        pt.show()


    #Plot: Pulsedpower VS Photocurrent
    def plot2(self):
        pt.plot(self.S_pulsedpower_integral,self.S_photocurrent_integral,'s')
        pt.xlabel('Laser Power Per Pulse [mJ]')
        pt.ylabel('Photocurrent Integral [nC] ')
        pt.title(self.name)
        pt.show()

    #Plot: CW power VS Photocurrent max & CW Power VS Pulsedpower max
    def plot3(self):
        fig, ax1 = pt.subplots()
        ax1.plot(self.S_cw_power,self.S_photocurrent_max_avg,'bs')
        ax2 = ax1.twinx()
        ax2.plot(self.S_cw_power,self.S_pulsedpower_max_avg,'rs')

        pt.show()


    #3D Scatter Plot: CW Power VS. Laser Power max VS. Photocurrent Max
    def plot4(self):
        fig = pt.figure()
        ax = fig.add_subplot(111,projection = '3d')
        
        x = self.S_cw_power
        y = self.S_pulsedpower_integral
        z = self.S_photocurrent_integral

        ax.scatter(x,y,z)
        ax.set_xlabel('cw power [\mu W]')
        ax.set_ylabel('Laser Power per Pulse [mJ]')
        ax.set_zlabel('Photocurrent Integral [nC]')
        pt.title(self.name)
        pt.grid()
        pt.show()

    #3D Surface Plot 
    def plot5(self):
        fig = pt.figure()
        ax = fig.gca(projection = '3d')
        
        x = self.S_cw_power
        y = self.S_pulsedpower_integral
        z = self.S_photocurrent_integral

        ax.plot_trisurf(x,y,z,cmap='cubehelix',linewidth=0.2)
        ax.set_xlabel('cw power [\mu W]')
        ax.set_ylabel('Laser Power per Pulse [mJ]')
        ax.set_zlabel('Photocurrent Integral [nC]')
        pt.title(self.name)

        pt.show()

    def plot6(self):
        fig = pt.figure()
        #plot
        #[cm.winter(i/10.) for i in range(10)]#[[0 84 159; 204 7 30; 87 171 39; 0 152 161; 97 33 88; 246 168 0]./255]
        #sns.palplot(sns.color_palette("hls", 8))
        #flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
        #sns.set_palette(flatui) #'hls'
        #figure(0) = pt.figure(dpi = 100)

        #Redefine the colors used to plot
        #mpl.rcParams['axes.color_cycle'] = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]#['#0000FF', '#00FF00', '#FF0000', '#00FFFF', 'FF00FF', 'FFFF00', '000000'] 
        
        for run in self.runs:
            x = run.pulsedpower_integral
            y = run.photocurrent_integral
            l = run.cw_power
            pt.plot(x,y,label='CW Power [\mu J]: '+str(l) )
            pt.xlabel('Laser Power per Pulse [mJ]')
            pt.ylabel('Photocurrent Integral [nC]')
        pt.legend()
        pt.title(self.name+': Photocurrent VS. Pulsedpower')
        #pt.rc('font',size=30)# not working if sns imported
        #pt.rc('font':'Helvetica Neue')
        #pt.rcParams.update({'font.size':24})
        pt.grid()
        pt.show()
    '''
    #something is not yet working because of the dimensions of the arrays
    def gaus_fits(self):
        fig, ax = pt.subplots(np.size(self.runs),sharex=True,sharey=True)
        for i in range(np.size(self.runs)):
            for j in self.runs[i].photocurrent_gausfit:
                x1 = self.runs[i].time
                y1 = self.runs[i].photocurrent_real
                x2 = self.runs[i].time
                y2 = gaus(self.runs[i].time,*j)
                print x1
                print y1
                print x2
                print y2
                ax[i].plot(x1,y1)
                ax[i].plot(x2,y2)
    '''


# This is a single run Object. It stores actual raw data
#from a single run with a certain CW Power
class single_run(sample):

    def __init__(self,filename):
        #call parent class and initialize data
        super(sample,self).__init__()

        self.filename = filename
        self.raw_data = np.loadtxt(filename)

        #data
        self.time             = None
        self.photocurrent    = None
        self.photodiode         = None
        self.cw_power        = None

        #converted to real units
        self.photocurrent_real = None
        self.pulsedpower_real = None


        # max and integral are linear so we can take the peak instead of the integral
        self.photocurrent_max = None
        self.pulsedpower_max = None

        #integral fits
        self.photocurrent_gausfit     = []
        self.pulsedpower_gausfit    = []

        #integrals
        self.photocurrent_integral     = None
        self.pulsedpower_integral     = None
        





    #this function simply plots all the raw data
    def plot_raw(self):
        pt.plot(np.transpose(self.raw_data))
        pt.show()




    #plots simple x and y of run #number
    def anyplot(self,x,y,xlabel,ylabel,title):


        #some color magic

        #plot
        #[cm.winter(i/10.) for i in range(10)]#[[0 84 159; 204 7 30; 87 171 39; 0 152 161; 97 33 88; 246 168 0]./255]
        #sns.palplot(sns.color_palette("hls", 8))
        #flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
        #sns.set_palette(flatui) #'hls'
        #figure(0) = pt.figure(dpi = 100)

        #Redefine the colors used to plot
        mpl.rcParams['axes.color_cycle'] = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]#['#0000FF', '#00FF00', '#FF0000', '#00FFFF', 'FF00FF', 'FFFF00', '000000'] 
        pt.figure(0)
        pt.plot(x,y,linewidth=2)
        pt.xlabel(xlabel)
        pt.ylabel(ylabel)
        pt.title(title +' '+ self.filename)
        pt.rc('font',size=30)# not working if sns imported
        #pt.rc('font':'Helvetica Neue')
        pt.rcParams.update({'font.size':24})
        pt.grid()
        pt.show()



    def plot_xy(self,y_value):
        #y_value can be: photocurrent, photodiode, photocurrent_real, pulsedpower_real

        #define x value
        x = np.transpose(self.time)*(10**9)
        xlabel         = 'time [ns]'


        #define y value
        if y_value == 'photocurrent':

            y = np.transpose(self.photocurrent)
            ylabel     = 'Photocurrent [V]'
            title     = 'Photocurrent: '

        elif y_value == 'photodiode':

            y = np.transpose(self.photodiode)
            ylabel     = 'Photodiode [V]'
            title     = 'Photodiode: '

        elif y_value == 'photocurrent_real':

            y = np.transpose(self.photocurrent_real*(10**6))
            ylabel     = 'Photocurrent[mA]'
            title     = 'Photocurrent: '

        elif y_value == 'pulsedpower_real':

            y = np.transpose(self.pulsedpower_real)
            ylabel     = 'Pulsed Laser Power[J/s]'
            title     = 'Pulsed Laser Power '


        else:
            print 'Error'

        # plot  data
        self.anyplot(x,y,xlabel,ylabel,title)


    # this function extracts time, photocurrent and photodiode data from the raw data
    def extract_data(self):
        #the data file format is

        #time
        #measurement t1 photocurrent
        #measurement t1 photodiode
        #measurement t2 photocurrent
        #measurement t2 photodiode
        #...

        self.time             = self.raw_data[0,:]          #first line
        self.photocurrent     = self.raw_data[1::2,:]    #every second line starting with the second line
        self.photodiode         = self.raw_data[2::2,:]     #every second line starting with the third line
        
        

    def subtract_noise(self,photocurrent_noise):
        
        self.photocurrent = self.photocurrent-photocurrent_noise
    
    
    #calculates the actual photocurrent in nanoAmperes and the  
    def return_real_units(self):

        #get the time in nanoseconds
        self.time = self.time * 10**9     #[ns]        
        
        '''
        photocurrent is measured in Ampere.
        I = U/(Gain * R_intern)
        '''
        #nano                     =10**9
        Gain                     = 5
        R_intern                 = 500 #Ohm
        self.photocurrent_real    = self.photocurrent/(Gain*R_intern) #[A]


        '''
        This is for the DET36 power meter
        pulsed laser power is measured in Watts
        V_out = P * R(lambda) * R_load
        P is laser power
        R(lambda) is return function and R(532)=0.199
        R_load is Resistance at the oscilloscope and 50 Ohms
        P = V_out/(R(lambda)*R_Load)*Gain
        [joule per second]
        '''
        
        Attenuator         = 10**-4    #[]
        Responsivity         = 0.199     #[A/W]
        R_load             = 50         #Ohm
        conversion_factor     = 1         #conversion factor from photodiode to energydetector 
        self.pulsedpower_real = self.photodiode/(Responsivity*R_load*Attenuator)*conversion_factor #J/s
        


    # fills the array with the maximum photocurrent / Pulsedpower
    def return_peaks(self):


        #remove offset
        pc_without_offset = self.return_without_offset(self.photocurrent_real,300)
        pp_without_offset = self.return_without_offset(self.pulsedpower_real,300)

        #get peaks
        self.photocurrent_max    = np.amax(pc_without_offset,axis=1)
        self.pulsedpower_max    = np.amax(pp_without_offset,axis=1)

        #add dimension back to array
        #self.photocurrent_max = np.expand_dims(self.photocurrent_max,axis=1)
        #self.pulsedpower_max = np.expand_dims(self.pulsedpower_max,axis=1)


        #print "time: ", np.shape(self.time)
        #print "photocurrent: ", np.shape(self.photocurrent)
        #print "photocurrent real: ", np.shape(self.photocurrent_real)
        #print "photocurrent without offset: ", np.shape(pc_without_offset)
        #use peakdetect        

    

    def return_integrals(self,gaus_guess,fit_only_a,fit_sigma):

        #Create empty arrays to store the gausfits and integrals
        self.photocurrent_integral = np.empty([np.size(self.photocurrent_real,0),1],dtype = float)
        self.pulsedpower_integral = np.empty([np.size(self.pulsedpower_real,0),1],dtype = float)

        if not(np.size(self.photocurrent_real,0) == np.size(self.pulsedpower_real,0)):
            print 'Photocurrent and Pulsedpower should have same length'

        #initial parameters for fitting function
        p_photocurrent     = gaus_guess[0]
        p_pulsedpower     = gaus_guess[1]

        #if we want to fit both a and sigma, we need to input both a and sigma
        if fit_sigma == True:
            def gaus_a(x,a,sigma):
                return gaus(x,a,p_photocurrent[1],sigma,p_photocurrent[3])
        #otherwise not
        else:
            def gaus_a(x,a):
                return gaus(x,a,p_photocurrent[1],p_photocurrent[2],p_photocurrent[3])
                
        for i in range(np.size(self.photodiode,0)):

            x = self.time
            #remember the time is in ns
            #nano  = 10**9
            #
            #milli = 10**6
                        
            
            y1 = self.photocurrent_real[i]
            y2 = self.pulsedpower_real[i]


            #if we only fit a, then we need to redefine the fit function and the initial guess
            if fit_only_a == True:
                
                fit_func = gaus_a
                #if we fit sigma too, we need to give 2 values as initial guesses
                if fit_sigma == True:
                    p_photocurrent_guess = [p_photocurrent[0],p_photocurrent[2]]
                #if we fit only a then 1 initial guess is good.
                else:
                    p_photocurrent_guess = [p_photocurrent[0]]
                
            else:
                fit_func = gaus
                p_photocurrent_guess = p_photocurrent
                
        
            #gaus(x,a,x0,sigma,y0):
            #fit photocurrent
            popt1_holder, pcov1_holder = curve_fit(fit_func,x,y1,p0=p_photocurrent_guess)#None,bounds = (0,[10.,1000.,500.]))
            #--> p[0]=a,p[1]=x0,p[2]=sigma,p[3]=y0
            
            #if we only fit a, then we need to rebuild the gaus parameters
            if fit_only_a == True:
                if fit_sigma == True:
                    popt1 = [popt1_holder[0],p_photocurrent[1],popt1_holder[1],p_photocurrent[3]]
                else:
                    popt1 = [popt1_holder[0],p_photocurrent[1],p_photocurrent[2],p_photocurrent[3]]
            else:
                popt1 = popt1_holder
            
            #Store the 
            self.photocurrent_gausfit.append(popt1)
            #save photocurrent integral
            #because time was mulitplied by 10**9 pu t into ns the integral is in nanoCoulomb
            self.photocurrent_integral[i] = gausint(popt1[0],popt1[2]) #[nanoCoulomb]

            popt2, pcov2 = curve_fit(gaus,x,y2,p0=p_pulsedpower)#None,bounds = (0,[10.,1000.,500.]))
            self.pulsedpower_gausfit.append(popt2)
            #becauset time was mulitplied by 10**9 and we want milliJoule
            self.pulsedpower_integral[i] = 10**-3*gausint(popt2[0],popt2[2]) #[mJ]            
            
            #print pcov1
            #print pcov2
            
            if False:
                if i%10==0:#self.photocurrent_integral[i]>0.0015:
                    #This plots the Raw Data and the Gaus Fits
                    #photocurrent            
                    pt.close()
                    fig = pt.figure()
                    ax = fig.add_subplot(111)
                    ym= gaus(self.time,popt1[0],popt1[1],popt1[2],popt1[3])
                    ax.plot(x,ym,c='r',label = 'best fit')
                    ax.scatter(x,y1,c='b',label='data photocurrent')
                    ax.legend()
                    fig.show()
                    raw_input('   ')
                    if self.photocurrent_integral[i]>0.0015:
                        print 'gaus_guess[0]: ', gaus_guess[0]
                        print 'p_photocurrent_guess: ',p_photocurrent_guess
                        print 'popt1: ',popt1
                        raw_input('press enter to continue...')
                    
            
            if False:
                pt.close()
                fig2 = pt.figure()
                ax = fig2.add_subplot(111)
                ym2= gaus(self.time,popt2[0],popt2[1],popt2[2],popt2[3])
                ax.plot(x,ym2,c='r',label = 'best fit')
                ax.scatter(x,y2,c='b',label='data pulsedpower')
                ax.legend()
                fig2.show()
            


            #pt.close("all")
            #input = raw_input('Enter: ')

            #if input == 'x':
            #    exit()
            
        

            




    def return_without_offset(self, data, left_boundary):

        #interesting. photocurrent integral is linear to max photocurrent
        #the same goes for photodiode integral and photodiode.
        #we can therfore just work with the peak value
        # this helps, because we do not need to define the left and right boundary for the integral.



        #get baseline
        offset_data = data[:,range(left_boundary)]
        offset         = np.average(offset_data,axis=1)
        offset         = np.expand_dims(offset,axis=1)
        #subtract baseline from data
        data         = data-offset
        return data
        #target = np.sum(data[:,range(left_boundary,right_boundary)],axis=1)


    #set CW power for specific run
    def set_cw(self,cw):

        self.cw_power = cw



#other functions

def gaus(x,a,x0,sigma,y0):
    return y0+a*np.exp(-(x-x0)**2/(2*sigma**2))

#int_{-inf,inf}[exp(-a*(x+b)^2)] = sqrt(PI/a)
#<->
def gausint(a,sigma):
    return  a * np.sqrt(np.pi*2*sigma**2)




#end of file data_analysis.py
