

library(celestial)
library(ProSpect)
library(stringr)
library(foreach)
library(magicaxis)



#Function to convert lookback time to redshift
lbt2z = function(z,lbt, ref){
  return(lbt - cosdistTravelTime(z = z, ref = ref))
  
} 

#skewed normal distribution but as a function of redshift
Evol_func = function(z, p){
  
  mSFR = p[1]
  mpeak = p[2]
  mperiod = p[3]
  mskew = p[4]
  
  snorm_norm_arg = ((z - mpeak) / mperiod) * exp(mskew)^asinh((z - mpeak) / mperiod) ##snorm 
  snorm  = exp(-1*snorm_norm_arg^2 / 2)
  
  mu = mpeak + abs(25 - mpeak)/2
  sigma = abs(25 - mpeak)/2
  
  trunc_part = 1 - 0.5*(1 + erf((z-mu)/(sigma*sqrt(2))))
  out = (snorm*trunc_part*mSFR)
  out[which(out < 0)] = 0
  return(out)
} 

#skewed normal function to generate CSFH
snorm=function(age,mSFR,mpeak,mperiod,mskew){
  mpeak=mpeak*1e9
  mperiod=mperiod*1e9
  age = age*1e9
  out = mSFR*dnorm(((age-mpeak)/mperiod)*(exp(mskew))^asinh((age-mpeak)/mperiod))*sqrt(2*pi)
  out[which(out < 0)] = 0
  return(out)
} 

#Covert the flux through a filter to an EBL prediction
ps2EBL_filt = function(f, freq, volweights){
  f = f * volweights[1:length(volweights)]
  sf = sum(f) # Sum epochs
  fout = f * 1e-26 # x1e-26 to get to W/m2/Hz (output is Jansky)
  fout = fout * 1e9 # to get to nW/m2/Hz
  fout = fout / (4*pi) # to get to nW/m2/Hz/sr
  
  
  return(sum(fout)*freq) 
} 
#Convert the sum of ProSpect SED's to an EBL model
ps2EBL = function(f, volweights, waveout){
  f = f * volweights[1:length(volweights)]
  sf = colSums(f) # Sum epochs
  fout = f * 1e-26 # x1e-26 to get to W/m2/Hz (output is Jansky)
  fout = fout * 1e9 # to get to nW/m2/Hz
  fout = fout / (4*pi) # to get to nW/m2/Hz/sr
  
  tempwave = 10^waveout
  tempfreq =  299792458/(tempwave) # Convert central wavelength from wavelength to freq
  sif= sf * tempfreq # Scale to nW/m2/sr
  return(fout) 
} 




ProSpectEBL = function(tau_birth = 1, 
                       tau_screen = 0.3, 
                       tau_AGN = 1,
                       pow_birth = -0.7,
                       pow_screen = -0.7,
                       pow_AGN = -0.7, 
                       alpha_SF_birth = 1,
                       alpha_SF_screen = 3, 
                       alpha_SF_AGN = 0, 
                       sparse = 5, 
                       speclib = NULL,
                       Dale = Dale_NormTot, 
                       AGN = Fritz, 
                       Dale_M2L_func = NULL,
                       returnall = TRUE, 
                       H0 = 67.8, 
                       OmegaM = 0.308,
                       OmegaL = 1 - OmegaM,
                       ref = '737', 
                       unimax = 1.38e+10, 
                       agemax = NULL,
                       LumDist_Mpc = NULL, 
                       addradio_SF = FALSE, 
                       addradio_AGN = FALSE, 
                       Te_SF = 10000,
                       ff_frac_SF = 0.1, 
                       ff_power_SF = -0.1, 
                       sy_power_SF = -0.8, 
                       Te_AGN = 10000,
                       ff_frac_AGN = 0.1, 
                       ff_power_AGN = -0.1, 
                       sy_power_AGN = -0.8, 
                       AGNct = 40, 
                       AGNrm = 60,
                       AGNan = 30, 
                       AGNta = 1, 
                       AGNal = 4, 
                       AGNbe = -0.5, 
                       AGNp = 1, 
                       AGNq = 1, 
                       Eb = 0, 
                       L0 = 2175.8,
                       LFWHM = 470, 
                       IGMabsorb = 0, 
                       waveout = seq(0,15,by = 0.01), 
                       #args unique to ProSpectEBL begin here.
                       agevec = seq(1e-5,13.5,0.06), 
                       zvec = NULL,
                       CSFH_parms = c(0.088,1.581,1.015,-0.299),
                       CSFH_vec = NULL,
                       AGN_parms = NULL,
                       AGN_vec = NULL,
                       type="redshift",
                       Zstart=1e-4,
                       Zfinal=0.02,
                       flts = NULL)
{
  
  
  #If age vector is provided without a corresponding redshift vector, backsolve to get
  #redshifts from ages. Computationally slow.
  if(is.null(zvec)){
    bins.redshift = foreach(lbt = agevec, .combine = "c") %do% newtonRaphson(f = lbt2z,
                                                                             x0 = 0,
                                                                             lbt = lbt,
                                                                             ref = ref)$root
    zvec = bins.redshift
  }
  
  else
    
    
    
    
  ch = cosdistTravelTime(zvec, ref = ref)
  agevec = cosdistTravelTime(zvec, ref=ref) * 1e9 # Create age vector
  age = agevec
  #if for some reason zvec (redshift vector) and age vector (agevec) are not the same length, stop.
  
  if(length(agevec) != length(zvec)){
    message("Redshift and age vectors are not equal in length! ")
    message("If provided, zvec must be the same length as agevec!")
    stop()
  }
  
  covolvec = cosdistCoVol(zvec, OmegaM = OmegaM, OmegaL=OmegaL, H0=H0, ref = ref) * 1e9 #  x1e9 to get to Mpc^3 (output is Gpc^3)
  volweights = (c(0,diff(covolvec)) + c(diff(covolvec),0))/2 # This evens out the volume weights around bin edges
  
  
  
  
  #Initialize Cosmology 
  tempsum = 0 #The variable we sum over
  
  #If CSFH is provided with snorm parameters and is a function of redshift, generate it. 
  if(is.null(CSFH_vec) & type == 'redshift'){
    SFRvec = Evol_func(zvec, p = CSFH_parms)
  }
  
  #Alternatively if CSFH is provided as a function of lookback time in years, generate CSFH.
  if(type == "lbt"){
    SFRvec = snorm(age = agevec/1e9, mSFR = CSFH_parms[1],
                   mpeak = CSFH_parms[2], mperiod = CSFH_parms[3],
                   mskew = CSFH_parms[4])
  }
  
  
  
  #Alternatively, if a constant star formation rate is provided, replicate it. 
  if(length(CSFH_vec) == 1){
    
    CSFH_vec = rep(CSFH_vec, length(agevec))
    
  }
  
  #Finally CSFH_vec is provided as a vector of values or a constant, overwrite parametric CSFH.
  if(!is.null(CSFH_vec)){
    SFRvec = CSFH_vec
  }
  
  tempSFH = approxfun(agevec, SFRvec, yleft = 0, yright = 0) # Generate our temporary SFH function to create our ZH
  
  
  #AGN parameters
  
  
  #if AGN lum is not provided or contains NULL, create vector of zeroes of the proper length to feed to ProSpectSED
  if(any(is.null(AGN_vec))){
    AGNlum = rep(0, length(agevec))
  }
  
  #If parameters for the CAGNH are provided, calculate AGNlum. Overwrites any constant/NULL AGNlum.
  if(any(!is.null(AGN_parms)) & type == 'redshift'){
    
    AGNlum = Evol_func(zvec, p = AGN_parms)
  }
  
  #If parameters for the CAGNH are provided, calculate AGNlum. Overwrites any constant/NULL AGNlum.
  if(any(!is.null(AGN_parms)) & type == 'lbt'){
    AGNlum = snorm(age = agevec/1e9, mSFR = AGN_parms[1],
                   mpeak = AGN_parms[2], mperiod = AGN_parms[3],
                   mskew = AGN_parms[4])
    
  }
  
  
  Zvec = Zfunc_massmap_box(agevec, massfunc=tempSFH, Zstart=Zstart, Zfinal=Zfinal) # Create out closed box ZH
  
  #Alternatively, if a constant AGN luminosity (of length 1) provided, replicate it
  if(length(AGN_vec) == 1){
    
    AGNlum = rep(AGN_vec, length(agevec))
    
  }
  
  
  #Check to make sure our vectors are all the same length.
  
  if(all(sapply(list(length(SFRvec),length(AGNlum),length(zvec),length(agevec)), function(x) x == length(agevec))) == FALSE){
    cat(paste(cbind('SFRvec',"AGNlum","zvec",'agevec')), "\n")
    cat(paste(cbind(length(SFRvec),length(AGNlum),length(zvec),length(agevec))), "\n")
    message("One or more vectors are not the same length. All vectors must be the same length.")
    stop()
  }
  
  
  #If provided, calculate the EBL as seen through filters provided in the flts argument.
  if(!is.null(flts)){
    
    freq_vals = c()  #Vector of pivot frequencies 
    filt_out = {} #list of approxfun() entries to pass to filtout in ProSpect
    pivot_wave = c() #get pivot wavelengths directly from getfilt()
    for(i in flts){filt_out=c(filt_out,list(approxfun(getfilt(i))))
    pivot_wave = append(pivot_wave, pivwavefunc(getfilt(i)))
    }
    freq_vals = 299792458/(pivot_wave/1e10)
    
    
    fluxz = foreach(i = 1:length(agevec), .combine='rbind')%do%{
      
      tempSFH = approxfun(agevec - agevec[i], SFRvec, yleft = 0, yright = 0) # Create a SFH from the current time to the end of the Universe
      tempZH = function(x,...){approxfun(agevec - agevec[i], Zvec, yleft = 0, yright = 0)(x)} # Create a ZH from the current time to the end of the Universe
      
      temp=ProSpectSED(SFH = SFHfunc, massfunc=tempSFH, z = zvec[i], agemax = max(agevec) - agevec[i], filters=NULL, # Run ProSpect SED for our target tempSFH and tempZH
                       alpha_SF_screen=alpha_SF_screen,alpha_SF_birth=alpha_SF_birth,
                       alpha_SF_AGN = alpha_SF_AGN, sparse = sparse,
                       Dale=Dale_NormTot, speclib=speclib, OmegaM = OmegaM, OmegaL=OmegaL, H0=H0, 
                       Z=tempZH, returnall=TRUE, addradio_SF=addradio_SF, ff_frac_SF = ff_frac_SF, 
                       Te_SF=Te_SF, ff_power_SF = ff_power_SF, sy_power_SF = sy_power_SF,
                       waveout = waveout, AGN = AGN, Te_AGN = Te_AGN, ff_frac_AGN = ff_frac_AGN,
                       ff_power_AGN = ff_power_AGN, sy_power_AGN = sy_power_AGN,
                       addradio_AGN = addradio_AGN, AGNlum = AGNlum[i], filtout = filt_out, AGNan = AGNan,
                       AGNal = AGNal, AGNta = AGNta, AGNct = AGNct, AGNbe = AGNbe, AGNrm = AGNrm, 
                       AGNp = AGNp, AGNq = AGNq, Eb = Eb, L0 = L0, LFWHM = LFWHM, IGMabsorb = IGMabsorb,  
                       pow_birth = pow_birth, pow_screen = pow_screen, pow_AGN = pow_AGN,
                       Dale_M2L_func = Dale_M2L_func, unimax = unimax)$Photom
      
      
    }
    
    
    output_data = matrix(data = 0, nrow = length(freq_vals),
                         ncol = 3)
    output_data = data.frame(output_data)
    colnames(output_data) = c("filter", "pivot_wave", "EBL_prediction")
    
    
    for(i in 1:length(freq_vals)){
      output_data$filter[i] = flts[i]
      output_data$pivot_wave[i] = pivot_wave[i]
      output_data$EBL_prediction[i] = ps2EBL_filt(fluxz[,i],299792458/(output_data$pivot_wave[i]/1e10),  volweights = volweights)
    }
    output_data$pivot_wave = output_data$pivot_wave / 1e4
    return(output_data)
  }
  
  #If no filter names are provided, generate an EBL model across all wavelengths.
  
  
  if(is.null(flts)){
    
    #now we will generate (and populate) empty matrices to populate with each of ProSpect's 
    #SED components (FinalFlux, StarsAtten, StarsUnatten, DustEmit, AGN)
    
    finalflux = matrix(0, nrow = length(agevec), ncol = length(waveout))
    stars_atten = matrix(0, nrow = length(agevec), ncol = length(waveout))
    StarsUnatten = matrix(0, nrow = length(agevec), ncol = length(waveout))
    DustEmit = matrix(0, nrow = length(agevec), ncol = length(waveout))
    AGN_flux = matrix(0, nrow = length(agevec), ncol = length(waveout))
    
    foreach(i = 1:length(agevec))%do%{
      
      tempSFH = approxfun(agevec - agevec[i], SFRvec, yleft = 0, yright = 0) # Create a SFH from the current time to the end of the Universe
      tempZH = function(x,...){approxfun(agevec - agevec[i], Zvec, yleft = 0, yright = 0)(x)} # Create a ZH from the current time to the end of the Universe
      
      
      temp=ProSpectSED(SFH = SFHfunc, massfunc=tempSFH, z = zvec[i], agemax = max(agevec) - agevec[i], filters=NULL, # Run ProSpect SED for our target tempSFH and tempZH
                       alpha_SF_screen=alpha_SF_screen,alpha_SF_birth=alpha_SF_birth,
                       alpha_SF_AGN = alpha_SF_AGN, sparse = sparse,
                       Dale=Dale_NormTot, speclib=speclib, OmegaM = OmegaM, OmegaL=OmegaL, H0=H0, 
                       Z=tempZH, returnall=TRUE, addradio_SF=addradio_SF, ff_frac_SF = ff_frac_SF, 
                       Te_SF=Te_SF, ff_power_SF = ff_power_SF, sy_power_SF = sy_power_SF,
                       waveout = waveout, AGN = AGN, Te_AGN = Te_AGN, ff_frac_AGN = ff_frac_AGN,
                       ff_power_AGN = ff_power_AGN, sy_power_AGN = sy_power_AGN,
                       addradio_AGN = addradio_AGN, AGNlum = AGNlum[i], filtout = NULL, AGNan = AGNan,
                       AGNal = AGNal, AGNta = AGNta, AGNct = AGNct, AGNbe = AGNbe, AGNrm = AGNrm, 
                       AGNp = AGNp, AGNq = AGNq, Eb = Eb, L0 = L0, LFWHM = LFWHM, IGMabsorb = IGMabsorb,  
                       pow_birth = pow_birth, pow_screen = pow_screen, pow_AGN = pow_AGN,
                       Dale_M2L_func = Dale_M2L_func, unimax = unimax)
      
      
      a = Lum2Flux(wave = temp$FinalLum$wave, lum = temp$FinalLum$lum, z = zvec[i], ref= ref)
      b = convert_wave2freq(flux_wave = a$flux, wave = a$wave)*1e23
      temp_fun = approxfun(x = a$wave, y = b)
      finalflux[i,] = temp_fun(10^waveout)
      #repeat this for all of the components...its ugly but functional
      # stars atten
      a = Lum2Flux(wave = temp$StarsAtten$wave, lum = temp$StarsAtten$lum, z = zvec[i], ref= ref)
      b = convert_wave2freq(flux_wave = a$flux, wave = a$wave)*1e23
      temp_fun = approxfun(x = a$wave, y = b)
      stars_atten[i,] = temp_fun(10^waveout)
      stars_atten[which(is.na(stars_atten))] = 0
      
      # stars unatten
      a = Lum2Flux(wave = temp$StarsUnAtten$wave, lum = temp$StarsUnAtten$lum, z = zvec[i], ref= ref)
      b = convert_wave2freq(flux_wave = a$flux, wave = a$wave)*1e23
      temp_fun = approxfun(x = a$wave, y = b)
      StarsUnatten[i,] = temp_fun(10^waveout)
      StarsUnatten[which(is.na(StarsUnatten))] = 0
      
      
      # Dust emission
      a = Lum2Flux(wave = temp$DustEmit$wave, lum = temp$DustEmit$lum, z = zvec[i], ref= ref)
      b = convert_wave2freq(flux_wave = a$flux, wave = a$wave)*1e23
      temp_fun = approxfun(x = a$wave, y = b)
      DustEmit[i,] = temp_fun(10^waveout)
      DustEmit[which(is.na(DustEmit))] = 0
      
      # AGN (if AGN lum is zero this is skipped,  for obvious reasons.)
      if(is.null(temp$AGN) == FALSE){
        a = Lum2Flux(wave = temp$AGN$wave, lum = temp$AGN$lum, z = zvec[i], ref= ref)
        b = convert_wave2freq(flux_wave = a$flux, wave = a$wave)*1e23
        temp_fun = approxfun(x = a$wave, y = b)
        AGN_flux[i,] = temp_fun(10^waveout)
        AGN_flux[which(is.na(AGN_flux))] = 0
      }
    }
    
    
    
    
    siflux = ps2EBL(finalflux, volweights = volweights, waveout = waveout)
    siflux_unatten = ps2EBL(StarsUnatten, volweights = volweights, waveout = waveout)
    siflux_atten = ps2EBL(stars_atten, volweights = volweights, waveout = waveout)
    siflux_dust = ps2EBL(DustEmit, volweights = volweights, waveout = waveout)
    siflux_AGN = ps2EBL(AGN_flux, volweights = volweights, waveout = waveout)
    
    
    tempwave = 10^waveout
    tempfreq =  299792458/(tempwave/1e10)
    tempsum = 0
    tempsum_unatten = 0
    tempsum_atten = 0
    tempsum_dust = 0
    tempsum_AGN = 0
    for(j in length(siflux[,1]):1){
      tempsum = tempsum + siflux[j,]
      tempsum_unatten = tempsum_unatten + siflux_unatten[j,]
      tempsum_atten = tempsum_atten + siflux_atten[j,]
      tempsum_dust = tempsum_dust + siflux_dust[j,]
      tempsum_AGN = tempsum_AGN + siflux_AGN[j,]
      
      
    }
    
    
    
    final_SED = cbind(tempwave/1e4, tempsum*tempfreq, tempsum_atten*tempfreq,
                      tempsum_unatten*tempfreq, tempsum_dust*tempfreq,
                      tempsum_AGN*tempfreq) #return wavelength in microns and EBL.
    final_SED = data.frame(final_SED)
    colnames(final_SED) = c("Wavelength", "Flux", "stars_atten",
                            'stars_unatten', 'dust', 'AGN')
    
    return(final_SED)
    
    
    
  }
  
  
}





