# ProSpectEBL_v0.5
A function that produces an EBL model from ProSpectSED.


ProSpectEBL Documentation 0.5

Dependencies - stringr - can be installed from the cran archive via

# install.packages("stringr")




Arguments:

type - String; Defines whether your cosmic star formation history and cosmic AGN activity history will be calculated 
as a function of redshift or lookback time in GYr. Must be either 'redshift' or 'lbt'. 

CSFH_parms - Numeric vector ; A vector of four pamameters that to use in the snorm(...) function that defines the cosmic star formation history in units of mSol yr^-1 Mpc^-3
as a function of lookback time or in Evol_func as function of redshift. In order the parameters are

mSFR - Numeric; The normalization for the underlying skewed normal distribution. 
	The highest star formation rate reached is equivalent to this value. Must, 
	of course, be positive and non-zero. Sensible values typically lie between 0.05-0.15 but any positive value is a valid input.

mpeak - Numeric; The look-back time (in GYr) or redshift at which the cosmic star formation rate peaks. Must be positive in either case and may
	return errors if mpeak is absurdly high (ie z=25 or lbt = 13.65). Goes without saying that if mpeak is given in lookback time it must be smaller 
	than the age of the universe. 
		
mperiod - Numeric; The parameter controlling the width of the underlying skewed normal distribution. Must be positive. Can be thought of as
	a standard deviation. A small value returns a narrow burst of star formation about mpeak, while larger values flatten out the SFH.
	
mskew - Numeric; The skew parameter of the underyling normal distribution. As a reminder, a value of zero makes the output equivalent
	to a standard normal distribution. Can, and often should be, positive. Negative values yield nonphysical star-formation histories which 
	skew towards high redshift (e.g. the CSFH is flipped and ends at zero while starting at high values.
	
See D'Silva et al. (2023) and Koushan et al. (2021) and references therein for discussion of the CSFH as a function of redshift or lookback time, respectively.

CSFH_vec - Numeric vector ; An option to pass in a vector as a cosmic star formation history in units mSol yr^-1 Mpc^-3. instead of calculating it via the snorm() functions.
Can be used to pass in non-parametric SFHs or SFHs calculated via a seperate user defined function. Default of NULL means that these values will be calculated
by the snorm function. If provided, overwrites the output of CSFH_parms. All values must be positive and this vector should be the same length as the agevec and/or zvec 
as each look-back time or redshift slice must have a corresponding star-formation rate. A single numeric value can also be provided and will be used at all time steps if for
any reason you want a constant star formation rate. Of course, all values must be positive.


AGN_parms - Numeric vector ; A vector of four parameters to represent the cosmic AGN luminosity history in units of Lsol erg s^-1. Uses the same 
skewed-normal distributions as CSFH_parms. See D'Silva et al. (2023) for a discussion of the CAGNLH. Goes without saying your AGN luminosity history 
parameters must be a function of the same variable as your star formation history, e.g. both must be functions of redshift or LBT. Mixing
the two will not return errors but will give nonsensical results. 
	

AGN_vec - Numeric vector ; An option to pass in a vector representing the AGN luminosity history in units of Lsol erg s^-1. Follows the same rules
as CSFH_vec and must be the same length as agevec, zvec, and CSFHvec if provided. As with CSFH_vec, if a single numeric value is provided it
will be taken as a constant AGN luminosity at all timesteps. Of course, all values must be positive.

Zstart - Numeric; The initial metallcity at which the model begins. Must be positive, see Zfunc_massmap_box for details.

Zfinal - Numeric; The final metallicity reached by the model at z=0. Must be positive, see Zfunc_massmap_box for details.

flts - String; A list of filters through which EBL predictions will be returned (e.g. 'g_VST', 'r_VST', 'i_VST'). 
Filters must exist in the ProSpect filter set through getfilt() or in EAZY filters. A default of NULL ignores filters and
returns a continuous EBL model as a function of wavelength.


agevec - Numeric vector; A vector of look-back times in GYr. Must be positive and all entries must be greater than zero,
sorted in ascending order from low to high LBT. An EBL prediction will be generated at each time-step and the final sum is returned.

zvec = Numeric vector; A vector of redshifts, must be in ascending order. If provided, must be the same length as agevec, otherwise the
redshifts corresponding to each LBT are calculated. If running ProSpectEBL repeatadly, calculating this redshift vector and passing
it in is recommended as back-solving for redshift from LBT is computationally slow.

*The following parameters are passed directly into ProSpectSED, which has already documented them in great detail. Additional notes,
if provided, explain details used by ProSpectEBL not included in the ProSpectSED documentation.


alpha_SF_screen

alpha_SF_birth

Dale

tau_screen

tau_birth

speclib

OmegaM

OmegaL

H0

returnall - Must be TRUE unlike ProSpectSED as ProSpectSED's components are used in all calculations.

addradio_SF

ff_frac_SF

Te_SF

waveout - Numeric; Desired output in log10 wavelength grid to use in Angstroms. Just a reminder since EBL studies often use units of microns but ProSpect uses Angstroms. 

addradio_AGN

AGNal

AGNta

AGNrm

AGNbe

AGNct

AGN - If NULL is provided, no AGN components will be calculted (see ProSpectSED). A value of NULL will mean AGN_parms and AGN_vec do nothing as
setting AGN to NULL disables AGN in ProSpectSED.

AGNan

ref 

ProSpectEBL returns a short list of output depending on whether a list of filters are provided to the *flts* argument. 

If a list of filters is not provided to *flts*, a table containing wavelength in microns versus EBL 
intensity (in nW m^-2 sr^-1) for each component of a ProSpectSED is returned.

If a list of filters to return predictions for is provided, instead of a table of wavelength versus EBL intensity, a table with three columns is returned.
The table contains the filter name (as defined in ProSpect), pivot wavelength (in microns) and EBL prediction (in nW m^-2 sr^-1). 



See Also

ProSpectSED, cosdistTravelTime, getfilt.



###Examples 

#Generate an EBL model which has no AGN component. CSFH from Bellstedt et al. (2020)
    
        fit_SFG = ProSpectEBL(CSFH_parms = c(0.069,10.94,2.074,0.394), Zfinal = 0.02, flts = NULL,
                          AGN_parms = NULL, type = 'lbt', agevec = seq(1e-2,13.5,0.06),
                          ref = '737')
     #Plot output. 

    magplot(fit_SFG$Wavelength, fit_SFG$Flux, log = "xy", xlim = c(0.1,1e5),
          type = 'l', side = c(1,2,3,4), ylim = c(1e-6,35), col = 'darkorchid1',xlab="Wavelength [microns]",
          ylab = expression("Intensity [nWm"^{-2} ~ "sr"^{-1} ~ "]"))
  
    lines(fit_SFG$Wavelength, fit_SFG$stars_atten, col = 'blue')
    lines(fit_SFG$Wavelength, fit_SFG$dust, col = 'chocolate4') 
    lines(fit_SFG$Wavelength, fit_SFG$AGN, col = 'red')
  
    legend(10^3.5, 100, legend = c("Total EBL", "Stars",
                                 'Dust', 'AGN'), bty = 'n', border = NULL, y.intersp = 0.35,
         text.col = c('darkorchid1', 'blue', 'chocolate4', 'red'), cex = 1.75)




#Generate an EBL model that has an AGN component. CSFH from Bellstedt et al. (2020) and CAGNH estimated from D'Silva et al 2023.
 
    fit_AGN = ProSpectEBL(CSFH_parms = c(0.069,10.94,2.074,0.394), Zfinal = 0.02, flts = NULL, AGN_vec = NULL,
                          AGN_parms = c(1e42,10.36546,1.75747,0.084666), type = 'lbt', agevec = seq(1e-2,13.5,0.06))

    #Plot Output
    
    magplot(fit_AGN$Wavelength, fit_AGN$Flux, log = "xy", xlim = c(0.1,1e5),
          type = 'l', side = c(1,2,3,4), ylim = c(1e-6,35), col = 'darkorchid1',xlab="Wavelength [microns]",
          ylab = expression("Intensity [nWm"^{-2} ~ "sr"^{-1} ~ "]"))
  
    lines(fit_AGN$Wavelength, fit_AGN$stars_atten, col = 'blue')
    lines(fit_AGN$Wavelength, fit_AGN$dust, col = 'chocolate4') 
    lines(fit_AGN$Wavelength, fit_AGN$AGN, col = 'red')
  
    legend(10^3.5, 100, legend = c("Total EBL", "Stars",
                                 'Dust', 'AGN'), bty = 'n', border = NULL, y.intersp = 0.35,
         text.col = c('darkorchid1', 'blue', 'chocolate4', 'red'), cex = 1.75)
  

#Generate EBL predictions though the VST/VISTA filter set. CSFH from Bellstedt et al. (2020). 	

	cols = c('darkorchid1', 'blue', 'green3', 
         'darkgoldenrod2', 'darkorange1', 'red', 'tan2','sienna4', 'black')
	fit_VST = try(ProSpectEBL(CSFH_parms = c(0.069,10.94,2.074,0.394), Zfinal = 0.02, 
                          flts = c('u_VST', 'g_VST', 'r_VST', 'i_VST', 'z_VST', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA'),
                          AGN_parms = NULL, type = 'lbt', agevec = seq(1e-2,13.5,0.06),
                          ref = '737'))


	magplot(fit_VST$pivot_wave, fit_VST$EBL_prediction, col = cols, log = 'xy',xlab="Wavelength [microns]",
        ylab = expression("Intensity [nWm"^{-2} ~ "sr"^{-1} ~ "]"), cex = 1.5, pch = 19,
        side = c(1,2,3,4))


	legend("bottomright", legend = c('u_VST', 'g_VST', 'r_VST', 'i_VST', 'z_VST', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA'),
        bty = 'n', border =NULL, text.col = cols)




                           
                           
