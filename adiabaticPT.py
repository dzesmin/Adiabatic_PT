#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import sys
import reader as rd

plt.ion()

# ======================================================================
# This code produces an adiabatic pressure and temperature profile  
# and generates 5 parameters that will be used as an initial guess
# for DEMC. The adiabatic profile is constrained based on the 
# planet effective temperature. The code reads a tep file, takes 
# the planet and stellar data, calculates planet's effective 
# temperature, produces an adiabatic profile (in the range
# of profiles in the literature) and plots it.
# The code takes 1 argument on the command line: "adibaticPT.py tepfile"
# Each run produces different profile so the user can pick the one
# that most suits his/her needs. Free parameters for the profile 
# are printed on the screen and returned by the function 'freeParams'.
# ======================================================================


# reads the tep file and calculates planet's effective temperature
def planet_Teff(tepfile):
     '''
     Calculates planetary effective temperature. Calls tep reader 
     and gets data needed for effective temperature calculation.
     The effective temperature is calculated assuming zero albedo, and 
     zero redistribution to the night side, i.e uniform dayside 
     redistribution.

     Parameters
     ----------
     tepfile: tep file, ASCII file
 
     Returns
     -------
     Teff: float 

     Example
     -------
     tepfile = "WASP-43b.tep"
     planet_Teff(tepfile)

     Revisions
     ---------
     2014-04-05 0.1  Jasmina Blecic, jasmina@physics.ucf.edu   Original version
     '''

     # opens tepfile to read and get data
     tep = rd.File(tepfile)

     # get stellar temperature
     stellarT = tep.getvalue('Ts')
     Tstar    = np.float(stellarT[0])

     # get stellar radius (in units of Rsun)
     stellarR = tep.getvalue('Rs')
     Rstar    = np.float(stellarR[0])

     # get semimajor axis in AU
     semimajor = tep.getvalue('a')
     a         = np.float(semimajor[0])

     # conversion to km
     AU   = 149597870.7   # km
     Rsun = 695500.0      # km

     # radius of the star and semimajor axis in km
     Rstar = Rstar * Rsun  # km
     a     = a * AU        # km

     # effective temperature of the planet 
     # Teff^4 = Teff*^4 * f * (Rstar/a)^2 * (1-A)
     # zero albedo, no energy redistribution to the night side A=0, f=1/2 
     Teff = Tstar * (Rstar/a)**0.5 * (1./2.)**0.25

     print '\n   The effective temperature of the planet is: ' + str(int(Teff)) + ' K'
     print

     return Teff


# generates free parameters
def freeParams(tepfile):
     '''
     This function produces free parameters for the adiabatic PT profile.
     It randomly samples empirically determined ranges for parameters.
     Temperature of the isothermal layer is constrained based on planet's
     effective temperature, assuming zero albedo and zero redistribution
     of the energy to the night side, and possible inclusion of the 
     spectral features in the range (1, 1.5)*Teff. 


     Parameters
     ----------
     tepfile: tep file, ASCII file
      
     Returns
     -------
     PT_params: 1D array of floats, free parameters

     Notes
     -----
     a1   , exponential factor
     a2   , exponential factor
     p1   , pressure at point 1
     p3   , pressure at point 3
     T3   , temperature of the layer 3, isothermal layer
     
     Revisions
     ---------
     2014-04-08 0.1  Jasmina Blecic, jasmina@physics.ucf.edu   Original version
     '''

     # takes Teff to constrain PT profile
     Teff = planet_Teff(tepfile)
     
     # sets empirical ranges of free parameters
     # temperature T3 constrained based on the planet dayside effective
     # temperature assuming zero albedo and zero redistribution to the night
     # side. The range set to Teff * (1, 1.5) to account for spectral features
     PTparams_range = np.array([
                            (0.99 , 0.999    ), # a1
                            (0.19 , 0.21     ), # a2
	                        (0.1  , 0.01     ), # p1
	                        (1    , 5        ), # p3 
	                        (Teff , Teff*1.5 ), # T3 
                                            ])

     # total number of free parameters
     noFreeParams = PTparams_range.shape[0]

     # calls random uniform for all free parameters
     base         = np.random.uniform(0., 1, noFreeParams)

     # free parameters of the PT profile
     PT_params    = base * (PTparams_range[:,1] - PTparams_range[:,0]) + PTparams_range[:,0]

     # prints in terminal free parameters
     print '\n   Free parameters of the adiabatic PT profile are: \n' + '        a1                a2               p1               p3             T3\n\n' + str(PT_params)
     print

     return PT_params


# calculates adiabatic PT profile
def PT_adiabatic(tepfile):
     '''
     This is a PT profile generator that uses similar methodology as derived in
     Madhusudhan and Seager 2009. It takes an equally spaced pressure array
     in log space from (100, 1e-5) bar and free parameters from freeParams function.
     It returns an adiabatic profile, a set of arrays for every layer in the
     atmosphere, a Gaussian smoothed temperature array, and a pressure.

     Parameters
     ----------
     tepfile: tep file, ASCII file
      
     Returns
     -------
     PT        : tuple of 1D arrays of floats
     T_smooth  : 1D array of floats
     p         : 1D array of floats

     Notes
     -----
     PT_NoInv      : tuple of temperatures and pressures for every layer in 
                     the atmosphere in non-inversion case
     T_smooth_Inv  : temperatures smoothed with Gaussian for inversion case
     p             : pressure in the atmosphere 

     Revisions
     ---------
     2014-04-05 0.1  Jasmina Blecic, jasmina@physics.ucf.edu   Original version
     '''

    # takes free parameters
     a1, a2, p1, p3, T3 = freeParams(tepfile)

     # sets pressures at the top and the bottom of the atmosphere
     bottom = 100    # bar
     top    = 1e-5   # bar
     p0     = top    # bar

     # set arbitrary number of levels in the atmosphere to 100
     noLevels = 100

     # the following set of equations derived using Equation 2
     # Madhusudhan and Seager 2009

     # temperature at point 1
     T1 = T3 - (np.log(p3/p1) / a2)**2

     # temperature at the top of the atmosphere
     T0 = T1 - (np.log(p1/p0) / a1)**2

     # error message when temperatures are < 0
     if T0<0 or T1<0 or T3<0:
          print 'T0, T1 and T3 temperatures are: ', T0, T1, T3
          raise ValueError('Input parameters give non-physical profile. Try again.')

     # log
     b  = np.log10(bottom)
     t  = np.log10(top)

     # equally spaced pressure space
     p = np.logspace(t, b, num=noLevels,  endpoint=True, base=10.0)

     # defining arrays for every part of the PT profile
     p_l1     = p[(np.where((p >= top) & (p < p1)))]
     p_l2_neg = p[(np.where((p >= p1)  & (p < p3)))]
     p_l3     = p[(np.where((p >= p3)  & (p <= bottom)))]

     # Layer 1 temperatures 
     T_l1 = (np.log(p_l1/p0) / a1)**2 + T0

     # Layer 2 temperatures decreasing part
     T_l2_neg = (np.log(p_l2_neg/p1) / a2)**2 + T1

     # Layer 3 temperatures
     T_l3 = np.linspace(T3, T3, len(p_l3))

     # concatenating all temperature arrays
     T_conc = np.concatenate((T_l1, T_l2_neg, T_l3))

     # PT profile
     PT = (T_l1, p_l1, T_l2_neg, p_l2_neg, T_l3, p_l3, T_conc, p, T0, T1, T3)

     # smoothing with Gaussian_filter1d
     sigma = 6
     T_smooth = gaussian_filter1d(T_conc, sigma, mode='nearest')

     return PT, T_smooth, p


# plots PT profiles
def plot_PT(tepfile):
     '''
     This function plots two figures:
     1.
     Adiabatic PT profile plotted part by part. It uses returned arrays from
     the PT_adiabatic function.
     2. 
     Smoothed PT profile without kinks on layer transitions.

     Parameters
     ----------
     tepfile: tep file, ASCII file
      
     Returns
     -------
     None

     Revisions
     ---------
     2014-04-08 0.1  Jasmina Blecic, jasmina@physics.ucf.edu   Original version
     '''  
    
     # generates adiabatic PT profile
     PT, T_smooth, p = PT_adiabatic(tepfile)

     # takes temperatures from PT generator
     T, T0, T1, T3 = PT[6], PT[8], PT[9], PT[10]

     # sets plots in the middle 
     minT= T0 * 0.9
     maxT= T3 * 1.1

     # plots raw PT profile
     plt.figure(1)
     plt.semilogy(PT[0], PT[1], '.', color = 'r'     )
     plt.semilogy(PT[2], PT[3], '.', color = 'b'     )
     plt.semilogy(PT[4], PT[5], '.', color = 'orange')
     plt.title('Adiabatic PT', fontsize=14)
     plt.xlabel('T [K]', fontsize=14)
     plt.ylabel('logP [bar]', fontsize=14)
     plt.xlim(minT  , maxT)
     plt.ylim(max(p), min(p))
     #plt.savefig('InitialAdiabaticPT.png', format='png')
     #plt.savefig('InitialAdiabaticPT.ps' , format='ps' )

     # plots Gaussian smoothing
     plt.figure(2)
     plt.semilogy(T_smooth, p, '-', color = 'b', linewidth=1)
     plt.title('Adiabatic PT Smoothed', fontsize=14)
     plt.xlabel('T [K]'     , fontsize=14)
     plt.ylabel('logP [bar]', fontsize=14)
     plt.ylim(max(p), min(p))
     plt.xlim(minT, maxT)
     #plt.savefig('InitialAdiabaticPTSmoothed.png', format='png')
     #plt.savefig('InitialAdiabaticPTSmoothed.ps' , format='ps' )

     return 


# executes the program
def main():
     '''
     This function sets the number of arguments that is called on the command
     line, and calls the PT generator function and the plotting function. 
     It keeps the plots on, until user closes them.

     Parameters
     ---------- 
     None
      
     Returns
     -------
     None

     Notes
     -----
     To call main, user needs to provide the following arguments:
     arg(1), tep file

     Revisions
     ---------
     2014-04-13  Jasmina Blecic, jasmina@physics.ucf.edu   Original version
     '''

     # counts number of arguments given
     noArguments = len(sys.argv)

     # prints usage if number of arguments different from 4
     if noArguments != 2:
        print '\nUsage: adiabaticPT.py <tepfile>'
        return
     
     # sets that the argument given is atmospheric file
     tepfile = sys.argv[1]

     print '\n\n    === Close the plots to end the program ===\n'

     # plots adiabatic PT profile and returns free parameters
     plot_PT(tepfile)

     # shows the plots until user closes them
     plt.show(block=True)


if __name__ == "__main__":
    main()


