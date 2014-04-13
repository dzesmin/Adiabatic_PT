

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

The code has the following functions:

planet_Teff(tepfile)
# reads the tep file and calculates planet's effective temperature

freeParams(tepfile)
# generates free parameters

PT_adiabatic(tepfile)
# calculates adiabatic PT profile

plot_PT(tepfile)
# plots PT profiles

main()
# executes the program

