import numpy as np
import multiprocessing as mp
from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc
#import exoplasim.surfacespecs
#import exoplasim.surfacespecs as spec
#import exoplasim.constants
#from exoplasim.constants import smws
import exoplasim.pyburn
import exoplasim.surfacespecs
import exoplasim.surfacespecs as spec
import exoplasim.constants
from exoplasim.constants import smws
from itertools import repeat
import exoplasim.colormatch


def _log(destination,string):
    if destination is None:
        if string=="\n":
            print()
        else:
            print(string)
    else:
        with open(destination,"a") as f:
            f.write(string+"\n")

def basicclouds(pressure,temperature,cloudwater):
    '''A basic cloud parameterization using T-dependent particle size distribution.
    
    This could be replaced with a different (better) cloud particle parameterization,
    but it should have the same call signature and return the same thing.
    
    Parameters
    ----------
    pressure : numpy.ndarray
        Pressure array for a column of the model
    temperature : numpy.ndarray
        Air temperatures for a column [K]
    cloudwater : numpy.ndarray
        Cloud water as a mass fraction [kg/kg] -- this is condensed water suspended in the cloud.
    
    Returns
    -------
    dict
        Dictionary of keyword arguments for setting clouds with an empirical particle size distribution.
    '''
    
    radius = {'H2O(c)':np.ones_like(pressure)*8e-4} #8 microns
    radius['H2O(c)'][temperature<=216.21] = (975.51*np.exp(0.05*(temperature[temperature<=216.21]-279.5))+0.66)*1e-4
    radius['H2O(c)'][temperature>216.21] = (175.44*np.exp(0.05*(temperature[temperature>216.21]-279.5))+34.45)*1e-4
    radius['H2O(c)'][temperature>240] = 20e-4 #20-micron water droplets
#Edwards et al 2007 doi:10.1016/j.atmosres.2006.03.002

    sigma_lnorm = 1.05
    
    return {'radius':radius, 'sigma_lnorm':sigma_lnorm}
    
def _adistance(lons,lats,tlon,tlat):
    '''Spherical law of Cosines to give angular distance in radians
    
    Parameters
    ----------
    lons : numpy.ndarray (1D) or float
        Longitudes in degrees against which to compare
    lats : numpy.ndarray (1D) or float
        Latitudes in degrees against which to compare
    tlon : float
        Target longitude in degrees
    tlat : float
        Target latitude in degrees
    
    Returns
    -------
    numpy.ndarray or float
        Angular distance to each in (lons,lats)
    '''
    
    rtlat = tlat*np.pi/180.
    rtlon = tlon*np.pi/180.
    rlons = lons*np.pi/180.
    rlats = lats*np.pi/180.
    distance = abs(np.arccos(np.sin(rtlat)*np.sin(rlats)+np.cos(rtlat)*np.cos(rlats)*np.cos(rlons-rtlon)))
    return distance

def _awidth(lons,lats,tlon,tlat):
    distances = np.sort(_adistance(lons,lats,tlon,tlat),axis=0)
    width = 0.5*(distances[1]+distances[2])
    return width
    
def _smooth(array,xcoords=None,xN=None,centerweight=0.95):
    '''Perform a conservative smoothing, such that the sum of either array or array*dx (if xcoords provided) is conserved.
    
    This routine keeps centerweight percent in each cell, and redistributes the rest into
    the neighboring cells.
    
    Parameters
    ----------
    array : 1D array
        Data array to be smoothed.
    xcoords : 1D array (optional)
        If provided, must have same shape as array. xN must also be provided.
    xN : float, optional
        Bookend value to use at the upper end of the array when computing dx.
    centerweight : float, optional
        The fraction of each cell that should be kept in the cell.
        
    Returns
    -------
    1D array
        Smoothed array
    '''
    padded = np.append(np.append([0.,],array),[0.,])
    edgeweight = 0.5*(1-centerweight)
    if xcoords is None:
        newarr = edgeweight*padded[:-2] + centerweight*padded[1:-1] + edgeweight*padded[2:]
    else:
        dx = np.diff(np.append(np.append(0.5*xcoords[0],0.5*(xcoords[:-1]+xcoords[1:])),xN))
        padded[1:-1] *= dx #brings us from density to mass
        newarr = (edgeweight*padded[:-2] + centerweight*padded[1:-1] + edgeweight*padded[2:])/dx
    return newarr

def _transitcolumn(atmosphere, pressures, psurf, temperature, humidity, clouds, 
                   gases_vmr, gascon, gravity, rplanet,
                   h2o_lines,cloudfunc,smooth,smoothweight,ozone,ozoneheight,ozonespread,num):
    '''Compute the transit spectrum through a single column.
    
    Parameters
    ----------
    atmosphere : Radtrans.Atmosphere
        Initialized Atmosphere object
    pressures : numpy.ndarray
        Pressure array in hPa, with the surface at the end of the array
    psurf : float
        Surface pressure in hPa
    temperature : numpy.ndarray
        Array of air temperatures, with the surface at the end of the array
    humidity : numpy.ndarray
        Specific humidity in the column [kg/kg]
    clouds : numpy.ndarray
        Cloud fraction in each cell
    gases_vmr : dict
        Dictionary of component gases and their volume mixing ratios
    gascon : float
        Specific gas constant
    gravity : float
        Surface gravity in SI units
    rplanet : float
        Planet radius in kilometers
    h2o_lines : {'HITEMP', 'EXOMOL'}
        Line list to use for H2O absorption. Options are 'HITEMP' or 'EXOMOL'.
    cloudfunc : function
        A routine which takes pressure, temperature, and cloud water content
        as arguments, and returns keyword arguments to be unpacked into calc_flux_transm.
    smooth : bool
        Whether or not to smooth humidity and cloud columns. As of Nov 12, 2021, it 
        is recommended that you use smooth=True for well-behaved spectra. This is a
        conservative smoothing operation, meaning the water and cloud column mass should
        be conserved--what this does is move some water from the water-rich layers into
        the layers directly above and below.
    smoothweight : float
        The fraction of the water in a layer that should be retained during smoothing.
        A higher value means the smoothing is less severe. 0.95 is probably the upper
        limit for well-behaved spectra.
    ozone : float
        Ozone column amount [cm-STP]
    ozoneheight : float
        Altitude of maximum ozone concentration [m]
    ozonespread : float
        Width of the ozone gaussian distribution [m]
    num : int, float
        Which number column this is (used for reporting)
        
        
    Returns
    -------
    numpy.ndarray
        Transmission radius in kilometers
    '''
    
    print("Doing column %d"%num)
    
    
    temperature = temperature[0]
    pressures = pressures[0]
    humidity = humidity[0]
    clouds = clouds[0]
    
    
    #Define vertical profiles with ghost stratosphere
    extp = np.append(np.geomspace(1.0e-6,0.5*pressures[0],num=20),pressures)
    try:
        extta = np.append(np.ones(20)*temperature[0],temperature)
    except:
        print(temperature.shape)
        raise
    exthus = np.append(np.zeros(20),humidity)
    extcl = np.append(np.zeros(20),clouds)
    
    if smooth:
        exthus = _smooth(exthus,xcoords=extp,xN=psurf,centerweight=smoothweight)
        extcl = _smooth(extcl,xcoords=extp,xN=psurf,centerweight=smoothweight)
    
    #Set up atmosphere
    atmosphere.setup_opa_structure(extp*1.0e-3)
    MMW = 8.31446261815324/gascon*1.0e3*np.ones_like(extp) #grams per mole
    
    #Set absorber fractions
    mass_fractions = {}
    for gas in gases_vmr:
        mmass = smws['m'+gas]
        mass_fractions[gas] = gases_vmr[gas] * mmass/MMW * (1-exthus-extcl)
        
    mass_fractions['H2O_'+h2o_lines] = exthus*(1-extcl)
    mass_fractions['H2O(c)'] = extcl
    kwargs = cloudfunc(extp,extta,extcl)
    
    #Compute ozone from parameterization
    mass_fractions['O3'] = np.zeros_like(extp)
    zfo3 = 100./2.14
    zo3t = ozone
    bo3 = ozoneheight
    co3 = ozonespread
    zh = 0.0
    zconst  = np.exp(-bo3/co3)
    ga = gravity*0.01 #SI units
    extdp = np.gradient(extp)*1e2
    for k in range(len(extp)-1,0,-1):
        zh -= extta[k]*gascon/ga*np.log(extp[k-1]/extp[k])
        zo3 = -(ozone + ozone*zconst)/(1.+np.exp((zh-bo3)/co3))+zo3t
        mass_fractions['O3'][k]=zo3*ga/(zfo3*extdp[k])
        zo3t -= zo3
    mass_fractions['O3'][0] = zo3t * ga / (zfo3 * extdp[0])
    
    if smooth:
        mass_fractions['O3'] = _smooth(mass_fractions['O3'],xcoords=extp,xN=psurf,centerweight=smoothweight)
    
    
    
    #Compute flux
    atmosphere.calc_transm(extta, mass_fractions, gravity, MMW, \
                           P0_bar=psurf*1.0e-3, R_pl=rplanet*1e5,**kwargs)
    
    return atmosphere.transm_rad*1e-5 #kilometers



def transit(output,transittimes,gases_vmr, gascon=287.0, gravity=9.80665, 
            rplanet=6.371e3,h2o_lines='HITEMP',num_cpus=4,cloudfunc=None,
            smooth=False,smoothweight=0.95,ozone=False,stepsperyear=11520.0,
            logfile=None):
    '''Compute transmission spectra for snapshot output
    
    This routine computes the transmission spectrum for each atmospheric column
    along the terminator, for each time in transittimes.
    
    Note: This routine does not currently include emission from atmospheric layers.
    
    Parameters
    ----------
    output : ExoPlaSim snapshot output
        Preferably opened with :py:func:`exoplasim.gcmt.load`.
    transittimes : list(int)
        List of time indices at which the transit should be computed.
    gases_vmr : dict
        Dictionary of gas species volume mixing ratios for the atmosphere
    gascon : float, optional
        Specific gas constant
    gravity : float, optional
        Surface gravity in SI units
    rplanet : float, optional
        Planet radius in km
    h2o_lines : {'HITEMP','EXOMOL'}, optional
        Either 'HITEMP' or 'EXOMOL'--the line list from which H2O absorption 
        should be sourced
    num_cpus : int, optional
        The number of CPUs to use
    cloudfunc : function, optional
        A routine which takes pressure, temperature, and cloud water content
        as arguments, and returns keyword arguments to be unpacked into calc_flux_transm.
    smooth : bool, optional
        Whether or not to smooth humidity and cloud columns. As of Nov 12, 2021, it 
        is recommended that you use smooth=True for well-behaved spectra. This is a
        conservative smoothing operation, meaning the water and cloud column mass should
        be conserved--what this does is move some water from the water-rich layers into
        the layers directly above and below.
    smoothweight : float, optional
        The fraction of the water in a layer that should be retained during smoothing.
        A higher value means the smoothing is less severe. 0.95 is probably the upper
        limit for well-behaved spectra.
        
    Returns
    -------
    Atmosphere,numpy.ndarray,numpy.ndarray,numpy.darray,numpy.ndarray
        pRT Atmosphere object, Wavelength in microns, array of all 
        transit columns with shape (ntimes,nterm,nfreq),
        where nterm is the number of terminator columns (time-varying), and nfreq
        is the number of frequencies in the spectrum, array of lon-lat coordinates for each
        transit specturm, array of spatial weights for each 
        column with the shape (ntimes,nterm) (for averaging), and the spatially-averaged
        transit spectrum, with shape (ntimes,nfreq). Transit radius is in km.
    '''
    _log(logfile,"===================================")
    _log(logfile,"| ExoPlaSim->petitRADTRANS Engine |")
    _log(logfile,"|   v1.0, Adiv Paradise (c) 2021  |")
    _log(logfile,"===================================")
    _log(logfile,"\n")
    _log(logfile,"--------Transit Spectrograph-------")
    _log(logfile,"\n")
    
    gravity *= 100.0 #SI->CGS
    
    _log(logfile,"Initializing petitRADTRANS atmosphere.....")
    
    atmosphere = Radtrans(line_species = ['H2O_'+h2o_lines,
                                          'CO2',
                                          'O3'],
                          rayleigh_species = ['N2', 'O2', 'CO2', 'H2', 'He'],
                          continuum_opacities = ['N2-N2', 'N2-O2','O2-O2',
                                                 'CO2-CO2','H2-H2','H2-He'],
                          wlen_bords_micron = [0.3, 20],
                          do_scat_emis = True,
                          cloud_species = ['H2O(c)_cd'],)
    
    atmosphere.hack_cloud_photospheric_tau = None
    
    if cloudfunc is None:
        cloudfunc = basicclouds
    
    transits = []
    weights = []
    terminatorlons = []
    terminatorlats = []
    meantransits = np.zeros((len(transittimes),len(atmosphere.freq)))
    
    _log(logfile,"Extracting model fields.....")
    
    lon = output.variables['lon'][:]
    lat = output.variables['lat'][:]
    lons,lats = np.meshgrid(lon,lat)
    nlon = len(lon)
    nlat = len(lat)
    lev = output.variables['lev'][:]
    tas = np.transpose(output.variables['ta'][:],axes=(0,2,3,1))
    h2os = np.transpose(output.variables['hus'][:],axes=(0,2,3,1))
    #clds = np.transpose(output.variables['cl'][:],axes=(0,2,3,1))
    dqls = np.transpose(output.variables['clw'][:],axes=(0,2,3,1))
    for idx,t in enumerate(transittimes):
        _log(logfile,"Configuring columns for timestamp %d corresponding to timestep %d....."%(t,output.variables['time'][t]))
        
        ts = output.variables['ts'][t,...].flatten()
        ps = output.variables['ps'][t,...].flatten()
        czen = output.variables['czen'][t,...]
        ta = np.reshape(tas[t,...],(nlat*nlon,len(lev)))
        hus = np.reshape(h2os[t,...],(nlat*nlon,len(lev)))
        #cld = np.reshape(clds[t,...],(nlat*nlon,len(lev)))
        dql = np.reshape(dqls[t,...],(nlat*nlon,len(lev)))
        pa = ps[:,np.newaxis]*lev[np.newaxis,:]
        
        #pa has shape [nlon*nlat,lev]
        
        #Extend down to the surface, using surface pressure, surface temperature, no clouds, and the same
        #humidity as the bottom vertical layer
        pa = np.concatenate([pa,ps[:,np.newaxis]],axis=1)
        hus = np.concatenate([hus,hus[:,-1][:,np.newaxis]],axis=1)
        ta = np.concatenate([ta,ts[:,np.newaxis]],axis=1)
        dql = np.concatenate([dql,np.zeros(len(dql))[:,np.newaxis]],axis=1)
        
        
        _log(logfile,"Identifying the terminator....")
        darkness = 1.0*(output.variables['czen'][t,...]==0.0)
        terminator_mask = darkness*(np.sqrt(np.gradient(darkness,axis=0)**2+\
                                            np.gradient(darkness,axis=1)**2)>0.0)
        terminator = np.argwhere(terminator_mask.flatten())
        psurf = ps[terminator]
        temp = ta[terminator,:]
        h2o = hus[terminator,:]
        clc = dql[terminator,:]
        press = pa[terminator,:]
        
        tlons = lons.flatten()[terminator]
        tlats = lats.flatten()[terminator]
        terminatorlons.append(tlons[:,0])
        terminatorlats.append(tlats[:,0])
        
        if ozone is False:
            a0o3 = 0.
            a1o3 = 0.
            aco3 = 0.
            bo3 = 20000.
            co3 = 5000.
            toffo3 = 0.0
        elif ozone is True:
            a0o3 = 0.25
            a1o3 = 0.11
            aco3 = 0.08
            bo3 = 20000.
            co3 = 5000.
            toffo3 = 0.25
        else:
            a0o3 = ozone["amount"]
            a1o3 = ozone["varlat"]
            aco3 = ozone["varseason"]
            toffo3 = ozone["seasonoffset"]
            bo3 = ozone["height"]
            co3 = ozone["spread"]
            
        dt = output.variables['time'][t] / stepsperyear
        rlats = tlats[:,0]*np.pi/180.
        o3 = a0o3+a1o3*abs(np.sin(rlats))+aco3*np.sin(rlats)*np.cos(2*np.pi*(dt-toffo3))
    
        
        widths = np.zeros_like(tlons[:,0])
        if num_cpus>1:
            args = zip(repeat(tlons),repeat(tlats),tlons,tlats)
            with mp.Pool(num_cpus) as pool:
                _w = pool.starmap(_awidth,args)
            for i,w in enumerate(_w):
                widths[i]=w
        else:
            for n in range(len(tlons)):
                widths[n] = _awidth(tlons,tlats,tlons[n],tlats[n])
        
        weights.append(widths)
        
        nterm = len(psurf)
        
        _log(logfile,"\n")
        _log(logfile,"Terminator mask applied; there are %d columns along the terminator."%nterm)
        
        transits.append(np.zeros((nterm,len(atmosphere.freq))))
        
        if num_cpus>1:
            _log(logfile,"\n")
            _log(logfile,"%d processes will be spun up; if this uses a substantial fraction\n"%num_cpus+
                         "of the computer's resources, it may become unresponsive for a while.")
            args = zip(repeat(atmosphere),press,psurf,temp,h2o,clc,repeat(gases_vmr),
                       repeat(gascon),repeat(gravity),repeat(rplanet),repeat(h2o_lines),
                       repeat(cloudfunc),repeat(smooth),repeat(smoothweight),o3,
                       repeat(bo3),repeat(co3),np.arange(nterm))
            with mp.Pool(num_cpus) as pool:
                spectra = pool.starmap(_transitcolumn,args)
            for i,column in enumerate(spectra):
                transits[-1][i,:] = column[:]
        else:
            _log(logfile,"\n")
            _log(logfile,"Running in single-process mode; this may take a while.")
            for i in range(nterm):
                print(i)
                transits[-1][i,:] = _transitcolumn(atmosphere,press[i,:],psurf[i],
                                                   temp[i,:],h2o[i,:],clc[i,:],
                                                   gases_vmr,gascon,gravity,
                                                   rplanet,h2o_lines,cloudfunc,
                                                   smooth,smoothweight,o3[i],bo3,co3,i)
        _log(logfile,"\n")
        _log(logfile,"All columns computed! Mean transit spectrum is now being computed.")
        _log(logfile,"\n")
        meantransits[idx,:] = np.average(transits[-1],axis=0,weights=widths)
    
    _log(logfile,"Repackaging into numpy arrays for export....")
    maxlen=0
    for transit in transits:
        maxlen=max(maxlen,len(transit[:,0]))
    nptransits = np.zeros((len(transittimes),maxlen,len(atmosphere.freq))) - 9999999.0 #finite fill value 
    npweights = np.zeros((len(transittimes),maxlen))
    npcoords = np.zeros((len(transittimes),maxlen,2)) + np.nan
    for n,transit in enumerate(transits):
        nptransits[n,:len(transit[:,0]),:] = transit[:,:]
        npweights[n,:len(transit[:,0])] = weights[n]
        npcoords[n,:len(transit[:,0]),0] = terminatorlons[n]
        npcoords[n,:len(transit[:,0]),1] = terminatorlats[n]
    return atmosphere,nc.c/atmosphere.freq*1e4,nptransits,npcoords,npweights,meantransits
        

def _imgcolumn(atmosphere, pressures, surface, temperature, humidity, clouds,
               gases_vmr, gascon, h2o_lines,
               gravity, Tstar, Rstar,
               starseparation,zenith,cloudfunc,smooth,smoothweight,ozone,
               ozoneheight,ozonespread):
    '''Compute the reflectance/emission spectrum for a column of atmosphere.
    
    Parameters
    ----------
    atmosphere : Radtrans.Atmosphere
        Initialized Atmosphere object
    pressures : numpy.ndarray
        Pressure array in hPa, with the surface at the end of the array
    surface : array-like
        Wavelength-array giving the surface reflectance, in decimal form.
    temperature : numpy.ndarray
        Array of air temperatures, with the surface at the end of the array
    humidity : numpy.ndarray
        Specific humidity in the column [kg/kg]
    clouds : numpy.ndarray
        Cloud fraction in each cell
    gases_vmr : dict
        Dictionary of component gases and their volume mixing ratios
    gascon : float
        Specific gas constant
    gravity : float
        Surface gravity in SI units
    h2o_lines : {'HITEMP', 'EXOMOL'}
        Line list to use for H2O absorption. Options are 'HITEMP' or 'EXOMOL'.
    Tstar : float
        Effective temperature of the parent star [K]
    Rstar : float
        Radius of the parent star in centimeters
    starseparation : float
        Distance between planet and star in centimeters
    zenith : float
        Angle between incident light and atmosphere normal vector in degrees.
    cloudfunc : function
        A routine which takes pressure, temperature, and cloud water content
        as arguments, and returns keyword arguments to be unpacked into calc_flux_transm.
    smooth : bool
        Whether or not to smooth humidity and cloud columns. As of Nov 12, 2021, it 
        is recommended that you use smooth=True for well-behaved spectra. This is a
        conservative smoothing operation, meaning the water and cloud column mass should
        be conserved--what this does is move some water from the water-rich layers into
        the layers directly above and below.
    smoothweight : float
        The fraction of the water in a layer that should be retained during smoothing.
        A higher value means the smoothing is less severe. 0.95 is probably the upper
        limit for well-behaved spectra.
    ozone : float
        Ozone column amount [cm-STP]
    ozoneheight : float
        Altitude of maximum ozone concentration [m]
    ozonespread : float
        Width of the ozone gaussian distribution [m]
    
    Returns
    -------
    numpy.ndarray
        Transmission radius in meters
    '''
    
    
    #Define vertical profiles with ghost stratosphere
    extp = np.append(np.geomspace(1.0e-6,0.5*pressures[0],num=20),pressures)
    extta = np.append(np.ones(20)*temperature[0],temperature)
    exthus = np.append(np.zeros(20),humidity)
    ##THIS BIT BC DRY LAYERS BREAK PRT
    #firstwater = np.argwhere(exthus>0.)[0][0]+1
    #exthus[firstwater:] = np.maximum(exthus[firstwater:],1.0e-5)
    ##REMOVE ONCE FIXED
    extcl = np.append(np.zeros(20),clouds)
    
    if smooth:
        exthus = _smooth(exthus,xcoords=extp,xN=psurf,centerweight=smoothweight)
        extcl = _smooth(extcl,xcoords=extp,xN=psurf,centerweight=smoothweight)
    
    #Set up atmosphere
    atmosphere.setup_opa_structure(extp*1.0e-3)
    MMW = 8.31446261815324/gascon*1.0e3*np.ones_like(extp) #grams per mole
    
    #Set absorber fractions
    mass_fractions = {}
    for gas in gases_vmr:
        mmass = smws['m'+gas]
        mass_fractions[gas] = gases_vmr[gas] * mmass/MMW * (1-exthus)
        
    mass_fractions['H2O_'+h2o_lines] = exthus
    kwargs = cloudfunc(extp,extta,extcl)
    
    #Compute ozone from parameterization
    mass_fractions['O3'] = np.zeros_like(extp)
    zfo3 = 100./2.14
    zo3t = ozone
    bo3 = ozoneheight
    co3 = ozonespread
    zh = 0.0
    zconst  = np.exp(-bo3/co3)
    ga = gravity*0.01 #SI units
    extdp = np.gradient(extp)*1e2
    for k in range(len(extp)-1,0,-1):
        zh -= extta[k]*gascon/ga*np.log(extp[k-1]/extp[k])
        zo3 = -(ozone + ozone*zconst)/(1.+np.exp((zh-bo3)/co3))+zo3t
        mass_fractions['O3'][k]=zo3*ga/(zfo3*extdp[k])
        zo3t -= zo3
    mass_fractions['O3'][0] = zo3t * ga / (zfo3 * extdp[0])
    
    if smooth:
        mass_fractions['O3'] = _smooth(mass_fractions['O3'],xcoords=extp,xN=psurf,centerweight=smoothweight)
    
    
    
    #Source hi-res surface spectrum from modelspecs, and make sure it matches the model
    #albedo
    
    atmosphere.reflectance = np.interp(nc.c/atmosphere.freq/1e-4,spec.wvl,surface)
    
    if zenith<=90.0: #We only need to compute the direct light contribution if the sun is up
        kwargs["theta_star":zenith]
        kwargs["Tstar":Tstar]
        kwargs["Rstar":Rstar]
        kwargs["semimajoraxis":starseparation]
    
    #Compute flux
    atmosphere.calc_flux(extta, mass_fractions, gravity, MMW,
                         geometry='non-isotropic',**kwargs)
    
    intensities = makeintensities(nc.c/atmosphere.freq/1e-4,atmosphere.flux*1e6)
    
    return atmosphere.flux*1e6,intensities #erg cm-2 s-1 Hz-1


def _lognorm(x):
    v = np.log10(np.maximum(x,1.0e-15))
    vmin = np.amin(v)
    v-=vmin
    vmax = np.amax(v)
    v/=vmax
    return v

def makeintensities(wvl,fluxes):
    '''Convert spectrum to (x,y,Y) intensities.
    
    Parameters
    ----------
    wvl : array-like
        Wavelengths in microns
    fluxes : array-like
        Spectrum in fluxes (units are arbitrary)
        
    Returns
    -------
    (float,float,float)
        (x,y,Y) tuple
    '''
    intensities = colormatch.makexyz(wvl*1.0e3,fluxes)
    return intensities
    
def makecolors(intensities):
    
    norms = _lognorm(intensities[...,2]/np.nanmax(intensities[...,2]))
    ogshape = intensities.shape
    flatnorms = norms.flatten()
    flatshape = len(flatnorms)
    flatintensities = np.reshape(intensities,(flatshape,3))
    colors = np.zeros((flatshape,3))
    for k in range(flatshape):
        colors[k,:] = colormatch.xyz2rgb(flatintensities[k,0],flatintensities[k,1],
                                         flatintensities[k,2]*flatnorms[k])
    colors /= np.nanmax(colors)
    colors = np.reshape(colors,ogshape)
    return colors
    
    
def image(output,imagetimes,gases_vmr, obsv_coords, gascon=287.0, gravity=9.80665, 
            Tstar=5778.0,Rstar=1.0,orbdistances=1.0,h2o_lines='HITEMP',
            num_cpus=4,cloudfunc=None,smooth=True,smoothweight=0.95,
            stellarspec=None,ozone=False,stepsperyear=11520.,logfile=None):
    '''Compute reflection+emission spectra for snapshot output
    
    This routine computes the reflection+emission spectrum for the planet at each
    indicated time.
    
    Note that deciding what the observer coordinates ought to be may not be a trivial operation.
    Simply setting them to always be the same is fine for a 1:1 synchronously-rotating planet,
    where the insolation pattern never changes. But for an Earth-like rotator, you will need to
    be mindful of rotation rate and the local time when snapshots are written. Perhaps you would
    like to see how things look as the local time changes, as a geosynchronous satellite might observe,
    or maybe you'd like to only observe in secondary eclipse or in quadrature, and so the observer-facing
    coordinates may not be the same each time.
    
    Parameters
    ----------
    output : ExoPlaSim snapshot output
        Preferably opened with :py:func:`exoplasim.gcmt.load`.
    imagetimes : list(int)
        List of time indices at which the image should be computed.
    gases_vmr : dict
        Dictionary of gas species volume mixing ratios for the atmosphere
    obsv_coords : list
        List of observer (lat,lon) coordinates for each
        observing time. Each item can either be a lat-lon tuple or a list of tuples
        (you can have multiple observers per snapshot). These are the surface coordinates 
        that are directly facing the observer. 
    gascon : float, optional
        Specific gas constant
    gravity : float, optional
        Surface gravity in SI units
    Tstar : float, optional
        Effective temperature of the parent star [K]
    Rstar : float, optional
        Radius of the parent star in solar radii
    orbdistances : float or numpy.ndarray, optional
        Distance between planet and star in AU
    h2o_lines : {'HITEMP','EXOMOL'}, optional
        Either 'HITEMP' or 'EXOMOL'--the line list from which H2O absorption 
        should be sourced
    num_cpus : int, optional
        The number of CPUs to use
    cloudfunc : function, optional
        A routine which takes pressure, temperature, and cloud water content
        as arguments, and returns keyword arguments to be unpacked into calc_flux_transm.
    smooth : bool, optional
        Whether or not to smooth humidity and cloud columns. As of Nov 12, 2021, it 
        is recommended that you use smooth=True for well-behaved spectra. This is a
        conservative smoothing operation, meaning the water and cloud column mass should
        be conserved--what this does is move some water from the water-rich layers into
        the layers directly above and below.
    smoothweight : float, optional
        The fraction of the water in a layer that should be retained during smoothing.
        A higher value means the smoothing is less severe. 0.95 is probably the upper
        limit for well-behaved spectra.
    stellarspec : array-like (optional)
        A stellar spectrum measured at the wavelengths in surfacespecs.wvl. If None, a
        blackbody will be used.
    ozone : bool or dict, optional
        True/False/dict. Whether or not forcing from stratospheric ozone should be included. If a dict
        is provided, it should contain the keys "height", "spread", "amount","varlat","varseason",
        and "seasonoffset", which correspond to the height in meters of peak O3 concentration, the 
        width of the gaussian distribution in meters, the baseline column amount of ozone in cm-STP, 
        the latitudinal amplitude, the magnitude of seasonal variation, and the time offset of the
        seasonal variation in fraction of a year. The three amounts are additive. To set a uniform, 
        unvarying O3  distribution, ,place all the ozone in "amount", and set "varlat" and 
        "varseason" to 0.
    stepsperyear : int or float, optional
        Number of timesteps per sidereal year. Only used for computing ozone seasonality.
        
    Returns
    -------
    Atmosphere, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray
        pRT Atmosphere object, wavelength in microns, Numpy array with dimensions (ntimes,nlat,nlon,nfreq),
        where ntimes is the number of output times, and nfreq
        is the number of frequencies in the spectrum, longitudes, latitudes, and an array with dimensions
        (ntimes,nfreq) corresponding to disk-averaged spectra, with individual contributions weighted by
        visibility and projected area.
    '''
    
    _log(logfile,"===================================")
    _log(logfile,"| ExoPlaSim->petitRADTRANS Engine |")
    _log(logfile,"|   v1.0, Adiv Paradise (c) 2021  |")
    _log(logfile,"===================================")
    _log(logfile,"\n")
    _log(logfile,"-------------Direct Imager----------")
    _log(logfile,"\n")
    
    #Compute zenith angle from czen
    #Figure out how to output spectra
    #Get surface type from model output
    #Turn off reflection for night-time grid cells
    
    gravity *= 100.0 #SI->CGS
    
    atmosphere = Radtrans(line_species = ['H2O_'+h2o_lines,
                                          'CO2',
                                          'O3'],
                          rayleigh_species = ['N2', 'O2', 'CO2', 'H2', 'He'],
                          continuum_opacities = ['N2-N2', 'N2-O2','O2-O2',
                                                 'CO2-CO2','H2-H2','H2-He'],
                          wlen_bords_micron = [0.3, 20],
                          do_scat_emis = True,
                          cloud_species = ['H2O(c)_cd'],)
    
    if cloudfunc is None:
        cloudfunc = basicclouds
    
    planckh =  6.62607015e-34
    boltzk = 1.380649e-23
    cc = 2.99792458e8
    
    if stellarspec is None:
        SIwvl = spec.wvl*1.0e-6 #um->m
        stellarspec = 1.0e-6 * 2*planckh*cc**2/SIwvl**5    \
                      /(np.exp(planckh*cc/(SIwvl*boltzk*Tstar))-1)
        # W sr-1 m-2 um-1
    lon = output.variables['lon'][:]
    lat = output.variables['lat'][:]
    lons,lats = np.meshgrid(lon,lat)
    nlon = len(lon)
    nlat = len(lat)
    
    images = np.zeros((len(imagetimes),nlat*nlon,len(atmosphere.freq)))
    photos = np.zeros((len(imagetimes),nlat*nlon,3))
    meanimages = np.zeros((len(imagetimes),len(atmosphere.freq)))
    
    lev = output.variables['lev'][:]
    tas = np.transpose(output.variables['ta'][:],axes=(0,2,3,1))
    h2os = np.transpose(output.variables['hus'][:],axes=(0,2,3,1))
    #clds = np.transpose(output.variables['cl'][:],axes=(0,2,3,1))
    dqls = np.transpose(output.variables['clw'][:],axes=(0,2,3,1))
    
    lsm = output.variables['lsm'][0,...].flatten()
    sea = 1.0-lsm #1 if sea, 0 if land
    
    ilons = lons.flatten()
    ilats = lats.flatten()
    
    lt1 = np.zeros(len(lat)+1)
    lt1[0] = 90
    for n in range(0,len(lat)-1):
        lt1[n+1] = 0.5*(lat[n]+lat[n+1])
    lt1[-1] = -90
    dln = np.diff(lon)[0]
    ln1 = np.zeros(len(lon)+1)
    ln1[0] = -dln
    for n in range(0,len(lon)-1):
        ln1[n+1] = 0.5*(lon[n]+lon[n+1])
    ln1[-1] = 360.0-dln
    
    lt1*=np.pi/180.0
    ln1*=np.pi/180.0
    
    darea = np.zeros((nlat,nlon))
    for jlat in range(0,nlat):
        for jlon in range(0,nlon):
            dln = ln1[jlon+1]-ln1[jlon]
            darea[jlat,jlon] = abs(np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*abs(dln)
    darea = darea.flatten()
    
    surfaces = [spec.modelspecs["groundblend"],
                spec.modelspecs["oceanblend"],
                spec.modelspecs["iceblend"]]
    
    sfcalbedo = surfaces[1][np.newaxis,:]*lsm[:,np.newaxis] + surfaces[0][np.newaxis,:]*(1-lsm)[:,np.newaxis]
    #sfcalbedo = sfcalbedo.flatten()
    
    for idx,t in enumerate(imagetimes):
        ts = output.variables['ts'][t,...].flatten()
        ps = output.variables['ps'][t,...].flatten()
        czen = output.variables['czen'][t,...]
        ta = np.reshape(tas[t,...],(nlat*nlon,len(lev)))
        hus = np.reshape(h2os[t,...],(nlat*nlon,len(lev)))
        #cld = np.reshape(clds[t,...],(nlat*nlon,len(lev)))
        dql = np.reshape(dqls[t,...],(nlat*nlon,len(lev)))
        pa = ps[:,np.newaxis]*lev[np.newaxis,:]
        
        #pa has shape [nlon*nlat,lev]
        
        #Extend down to the surface, using surface pressure, surface temperature, no clouds, and the same
        #humidity as the bottom vertical layer
        pa = np.concatenate([pa,ps[:,np.newaxis]],axis=1)
        hus = np.concatenate([hus,hus[:,-1][:,np.newaxis]],axis=1)
        ta = np.concatenate([ta,ts[:,np.newaxis]],axis=1)
        dql = np.concatenate([dql,np.zeros(len(dql))[:,np.newaxis]],axis=1)
        
        if ozone is False:
            a0o3 = 0.
            a1o3 = 0.
            aco3 = 0.
            bo3 = 20000.
            co3 = 5000.
            toffo3 = 0.0
        elif ozone is True:
            a0o3 = 0.25
            a1o3 = 0.11
            aco3 = 0.08
            bo3 = 20000.
            co3 = 5000.
            toffo3 = 0.25
        else:
            a0o3 = ozone["amount"]
            a1o3 = ozone["varlat"]
            aco3 = ozone["varseason"]
            toffo3 = ozone["seasonoffset"]
            bo3 = ozone["height"]
            co3 = ozone["spread"]
            
        dt = output.variables['time'][t] / stepsperyear
        rlats = ilats*np.pi/180.
        o3 = a0o3+a1o3*abs(np.sin(rlats))+aco3*np.sin(rlats)*np.cos(2*np.pi*(dt-toffo3))
    
        
        zenith = np.arccos(czen)*180./np.pi
        darkness = 1.0*(czen==0.0)
        nightside = darkness*(np.sqrt(np.gradient(darkness,axis=0)**2+\
                                      np.gradient(darkness,axis=1)**2)==0.0)
        zenith[nightside>0.5] = 91.0
        zenith = zenith.flatten()
        
        albedo = output.variables['alb'][t,...].flatten()
        ice = output.variables['sit'][t,...]+output.variables['snd'][t,...]
        ice = ice.flatten()
        #icemap = 2.0*(ice>0.001) #1 mm probably not enough to make everything white
        ice = np.minimum(ice/0.02,1.0) #0-1 with a cap at 2 cm of snow
        #sfctype = np.maximum(sea,icemap).astype(int)
           # Sea with no ice:  1
           # Sea with ice:     2
           # Land with no ice: 0   (since both are 0)
           # Land with ice:    2
        #surfaces = []
        #for n in range(len(sfctype)):
        #    surfaces.append({'type':sfctype[n],'albedo':albedo[n]})
            
        surfspecs = (sfcalbedo*(1-output.variables['sic'][t,...].flatten())[:,np.newaxis]+
                     surfaces[2][np.newaxis,:]*(output.variables['sic'][t,...].flatten())[:,np.newaxis])
        surfspecs = (surfspecs*(1-ice)[:,np.newaxis] + surfaces[2][np.newaxis,:]*ice[:,np.newaxis])
        
        for n in range(len(albedo)):
            bol_alb = np.trapz(stellarspec*surfspecs[n,:],x=spec.wvl)/ \
                      np.trapz(stellarspec,x=spec.wvl)
            fudge_factor = albedo[n]/bol_alb
            surfspecs[n,:] *= fudge_factor
            
        viewangles = []
        try:
            test = len(obsv_coords[idx][0]) #If this fails, then we only have one set of coordinates
            
            for coord in obsv_coords[idx]:
                view = _adistance(ilons,ilats,coord[1],coord[0])
                view[view>=np.pi/2.] = np.pi/2.
                viewangles.append(view)
        except:
            view = _adistance(ilons,ilats,obsv_coords[idx][1],obsv_coords[idx][0])
            view[view>=np.pi/2.] = np.pi/2.
            viewangles.append(view)
        viewangles = _adistance(ilons,ilats,obsv_coords[idx][1],obsv_coords[idx][0])
        
        projectedareas = []
        for view in viewangles:
            projectedareas.append(np.cos(view)*darea)
        
        if num_cpus>1:
            args = zip(repeat(atmosphere),pa,surfspecs,ta,hus,dql,repeat(gases_vmr),
                       repeat(gascon),repeat(h2o_lines),repeat(gravity),
                       repeat(Tstar),repeat(Rstar*nc.r_sun),repeat(starseparation*nc.AU),
                       zenith,repeat(cloudfunc),repeat(smooth),repeat(smoothweight),o3,
                       repeat(bo3),repeat(co3))
            with mp.Pool(num_cpus) as pool:
                spectra = pool.starmap(_imgcolumn,args)
            for i,column in enumerate(spectra):
                images[idx,i,:] = column[0][:]
                photos[idx,i,:] = column[1][:]
        else:
            for i in range(nterm):
                column = _imgcolumn(atmosphere,pa[i,:],surfspecs[i,:],
                                    ta[i,:],hus[i,:],dql[i,:],
                                    gases_vmr,gascon,h2o_lines,
                                    gravity,Tstar,Rstar*nc.r_sun,
                                    starseparation*nc.AU,zenith[i],
                                    cloudfunc,smooth,smoothweight,o3[i],
                                    bo3,co3)
                images[idx,i,:] = column[0][:]
                photos[idx,i,:] = column[1][:]
        for idv,view in enumerate(projectedareas):
            meanimages[idx,idv,:] = np.average(images[idx,...],axis=0,weights=view)
    
    images = np.reshape(images,(len(imagetimes),nlat,nlon,len(atmosphere.freq)))
    photos = np.reshape(photos,(len(imagetimes),nlat,nlon,3))
    return atmosphere,nc.c/atmosphere.freq*1e4,images,photos,lon,lat,meanimages

_postmetadata = {"images":["image_spectra_map","erg cm-2 s-1 Hz-1"],
                 "spectra":[["image_spectrum_avg","erg cm-2 s-1 Hz-1"],
                            ["transit_spectrum_avg","km"]],
                 "colors":["true_color_image","RGB"],
                 "transits":["transit_spectra_map","km"],
                 "weights":["transit_column_width","degrees"],
                 "lat":["latitude","degrees"],
                 "lon":["longitude","degrees"],
                 "wvl":["wavelength","microns"],
                 "time":["obsv_time_index","indices"]}

def _writecsvs(filename,variables,extension=None,logfile=None):
    '''Write CSV output files
    
    Files are placed in a subdirectory named from the filename naming pattern (stripping off the extension).
    
    Parameters
    ----------
    filename : str
        Filename pattern
    variables : dict
        Dictionary of variable data arrays
    meta : dict
        Dictionary of metadata fields for associated variables
    extension : str, optional
        File extension to use for individual files
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
        
    Returns
    -------
    list, str 
        List of paths to output files, and the containing directory.
    '''
    idx = filename[::-1].find(".")+1 #Index of last period separator (negative)
    dirname = filename[:-idx]
    if dirname[-4:]==".tar":
        dirname = dirname[:-4]
        idx+=4
    if extension is None:
        extension = filename[-idx:]
    fname = dirname
    if "/" in dirname:
        fname = dirname.split("/")[-1]
    os.system("mkdir %s"%dirname) #Create a subdirectory that just omits the file extension
    files = []
    if "images" in variables:
        dimvars = ["lat","lon","wvl","time"]
    else:
        dimvars = ["wvl","time"]
    keyvars = list(variables.keys())
    for var in dimvars:
        outname =  "%s/%s_%s%s"%(dirname,fname,var,extension)
        np.savetxt(outname,variables[var].astype("float32"),
                   header=(str(len(variables[var]))+",|||,"
                          +','.join(_postmetadata[var])),delimiter=',')
        keyvars.remove(var)
        files.append(outname)
    maxlen = 0
    for var in keyvars:
        maxlen = max(maxlen,len("%s/%s_%s%s"%(dirname,fname,var,extension)))
    maxlen+=1
    for var in keyvars:
        #This creates e.g. most_output/most_output_ts.csv if filename was most_output.csv
        shape = variables[var].shape
        dim2 = shape[-1]
        dim1 = int(np.prod(shape[:-1]))
        var2d = np.reshape(variables[var],(dim1,dim2)) #np.savetxt can only handle 2 dimensions
        outname =  "%s/%s_%s%s"%(dirname,fname,var,extension)
        if var=="spectra":
            if "images" in keyvars:
                meta = _postmetadata[var][0]
            else:
                meta = _postmetadata[var][1]
        else:
            meta = _postmetadata[var]
        try:
            np.savetxt(outname,var2d.astype("float64"),
                   header=(','.join(np.array(shape).astype(str))+",|||,"
                          +','.join(meta)),delimiter=',')
        except:
            print(",".join(np.array(shape).astype(str)))
            print(",|||,")
            print(meta)
            print(','.join(meta))
            print(','.join(np.array(shape).astype(str))+",|||,"
                          +','.join(meta))
            raise
        #The original shape of the array to which it should be reshaped on unpacking is in the header,
        #with the actual metadata separated from the shape by '|||'
        files.append(outname)
        try:
            writeline = "Writing %8s to %"+str(maxlen)+"s\t....... %d timestamps"
            _log(logfile,writeline%(var,outname,variables[var].shape[0]))
        except:
            _log(logfile,"Writing %8s to %s"%(var,filename))
    return files,dirname+"/"   
    
def _csv(dataset,filename="most_output.tar.gz",logfile=None,extracompression=False):
    '''Write a dataset to CSV/TXT-type output, optionally compressed.
    
    If a tarball format (e.g. \*.tar or \*.tar.gz) is used, output files will be packed into a tarball.
    gzip (.gz), bzip2 (.bz2), and lzma (.xz) compression types are supported. If a tarball format is 
    not used, then accepted file extensions are .csv, .txt, or .gz. All three will produce a directory
    named following the filename pattern, with one file per variable in the directory. If the .gz extension
    is used, NumPy will compress each output file using gzip compression. 
    
    Files will only contain 2D
    variable information, so the first N-1 dimensions will be flattened. The original variable shape is
    included in the file header (prepended with a # character) as the first items in a comma-separated
    list, with the first non-dimension item given as the '|||' placeholder. On reading variables from these
    files, they should be reshaped according to these dimensions. This is true even in tarballs (which 
    contain CSV files).
    
    Parameters
    ----------
    dataset : dict
        A dictionary of outputs as generated from :py:func:`pyburn.dataset()<exoplasim.pyburn.dataset>`
    filename : str, optional
        Path to the output file that should be written. This will be parsed to determine output type.
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
    extracompression : bool, optional
        If True, then component files in tarball outputs will be compressed individually with gzip, 
        instead of being plain-text CSV files.
        
    Returns
    -------
    tuple or str
        If non-tarball output was used, a tuple containing a list of paths to output files, and a string
        giving the name of the output directory. If tarball output was used, a relative path to the tarball.
    '''
    
    variables = {}
    
    fileparts = filename.split('.')

    if fileparts[-2]=="tar": #We will be doing compression
        import tarfile
        ext = ".csv"
        if extracompression:
            ext = ".gz"
        files,dirname = _writecsvs(filename,dataset,extension=ext,logfile=logfile)
        namelen = len(filename)
        with tarfile.open(filename,"w:%s"%fileparts[-1]) as tarball:
            
            maxlen = 0
            for tfile in files:
                maxlen = max(maxlen,len(tfile))
            for var in files:
                varname = var
                if len(varname.split("/"))>2:
                    varname = "/".join(varname.split("/")[-2:])
                tarball.add(var,arcname=varname)
                writeline = "Packing %"+str(maxlen)+"s in %"+str(namelen)+"s"
                try:
                    _log(logfile,(writeline+" ....... %d timestamps")%(var,filename,
                                                                       dataset[var].shape[0]))
                except:
                    _log(logfile,writeline%(var,filename))
        os.system("rm -rf %s"%dirname)
        return filename
        
    elif fileparts[-1]=="tar": #We're still making a tarball, but it won't be compressed
        import tarfile
        ext = ".csv"
        if extracompression:
            ext = ".gz"
        files,dirname = _writecsvs(filename,variables,extension=ext,logfile=logfile)
        namelen = len(filename)
        
        with tarfile.open(filename,"w") as tarball:
            
            maxlen = 0
            for tfile in files:
                maxlen = max(maxlen,len(tfile))
            for var in files:
                tarball.add(var)
                writeline = "Packing %"+str(maxlen)+"s in %"+str(namelen)+"s"
                try:
                    _log(logfile,(writeline+"\t....... %d timestamps")%(var,filename,
                                                                   dataset[var].shape[0]))
                except:
                    _log(logfile,writeline%(var,filename))
        os.system("rm -rf %s"%dirname)
        return filename
        
    else: #Just a collection of CSV/TXT-type files in a subdirectory, which may be individually-compressed.
        #These files can have .txt, .csv, or .gz file extensions.
        files,dirname = _writecsvs(filename,variables,logfile=logfile)
        return files,dirname
    
def _npsavez(filename,dataset,logfile=None):

    meta = {}
    for var in _postmetadata:
        if var=="spectra" and "images" in dataset:
            meta[var] = _postmetadata[var][0]
        elif var=="spectra" and "transits" in dataset:
            meta[var] = _postmetadata[var][1]
        else:
            meta[var] = _postmetadata[var]
    metafilename = filename[:-4]+"_metadata.npz"
    
    np.savez_compressed(metafilename,**meta)
    np.savez_compressed(filename,**dataset)
    
    return metafilename,filename

def _netcdf(filename,dataset,logfile=None):
    
    import netCDF4 as nc
    
    _log(logfile,"Writing petitRADTRANS output to %s."%filename)
    _log(logfile,"Initializing file and dimensions....")
    
    latitude  = dataset["lat"]
    longitude = dataset["lon"]
    time      = dataset["time"]
    wvls      = dataset["wvl"]
    
    ntimes = len(time)
    nwvls  = len(wvls)
    
    ncd = nc.Dataset(filename, 'w', format="NETCDF4")
    
    t0 = 0
    t1 = ntimes
    
    twoD = True
    
    if "images" in dataset:
        twoD = True
        _log(logfile,"Direct imaging data detected; writing output with 2D maps")
    else: #Transit spectra
        twoD = False
        _log(logfile,"Transit data detected; writing output with 1D terminator maps")
        
    dt = ncd.createDimension("time", ntimes)
    times = ncd.createVariable("time","f4",("time",), zlib=True, least_significant_digit=6)
    wvl = ncd.createDimension("wvl", nwvls)
    wavelengths = ncd.createVariable("wvl","f4",("wvl",), zlib=True, least_significant_digit=6)

    if twoD:
        nlats  = len(latitude)
        nlons  = len(longitude)
        lat = ncd.createDimension("lat", nlats)
        lon = ncd.createDimension("lon", nlons)
        latitudes  = ncd.createVariable("lat","f4",("lat",), zlib=True, 
                                        least_significant_digit=6)
        longitudes = ncd.createVariable("lon","f4",("lon",), zlib=True, 
                                        least_significant_digit=6)
    else:
        ncoords = latitude.shape[1]
        coord = ncd.createDimension("coord", ncoords)
        latitudes  = ncd.createVariable("lat","f4",("time","coord"), zlib=True, 
                                        least_significant_digit=6)
        longitudes = ncd.createVariable("lon","f4",("time","coord"), zlib=True, 
                                        least_significant_digit=6)
    
    ncd.set_auto_mask(False)
    times.set_auto_mask(False)
    latitudes.set_auto_mask(False)
    longitudes.set_auto_mask(False)
    wavelengths.set_auto_mask(False)
    
    times.units = "indices"
    latitudes.units = "degrees"
    longitudes.units= "degrees"
    wavelengths.units="microns"
    
    times[:]  = time[:].astype("float32")
    latitudes[:] = latitude[:].astype("float32")
    longitudes[:] = longitude[:].astype("float32")
    wavelengths[:] = wvls[:].astype("float32")
    
    times.axis = 'T'
    latitudes.axis = 'Y'
    longitudes.axis = 'X'
    wavelengths.axis= 'Z'
    
    times.standard_name = "obsv_time_index"
    latitudes.standard_name = "latitude"
    longitudes.standard_name = "longitude"
    wavelengths.standard_name= "wavelength"
    
    times.long_name = "obsv_time_index"
    latitudes.long_name = "latitude"
    longitudes.long_name = "longitude"
    wavelengths.long_name= "wavelength"
    
    if twoD:
        _log(logfile,"Writing 2D image data for each observation.....")
        datavar = dataset["images"]
        dvarmask = datavar[np.isfinite(datavar)]
        dvarmask = dvarmask[dvarmask!=0]
        if len(dvarmask)==0:
            lsd = None
        else:
            lsd = max(int(round(abs(np.log10(abs(dvarmask).min()))+0.5))+6,6) #get decimal place of smallest value
            if abs(dvarmask).min()>=1.0:
                lsd=6
        images = ncd.createVariable("images","f8",("time","lat","lon","wvl"),
                                    zlib=True,least_significant_digit=lsd)
        images.set_auto_mask(False)
        images[:] = datavar[:]
        images.units = "erg cm-2 s-1 Hz-1"
        images.standard_name = "image_spectra_map"
        images.long_name = "image_spectra_map"
        
        _log(logfile,"Writing 2D true-color data for each observation.....")
        datavar = dataset["colors"]
        dvarmask = datavar[np.isfinite(datavar)]
        dvarmask = dvarmask[dvarmask!=0]
        photos = ncd.createVariable("colors","f4",("time","lat","lon","RGB"),
                                    zlib=True,least_significant_digit=6)
        photos.set_auto_mask(False)
        photos[:] = datavar[:]
        photos.units = "RGB"
        photos.standard_name = "true_color_image"
        photos.long_name = "true_color_image"
        
        _log(logfile,"Writing disk-averaged image data for each observation.....")
        datavar = dataset["spectra"]
        dvarmask = datavar[np.isfinite(datavar)]
        dvarmask = dvarmask[dvarmask!=0]
        if len(dvarmask)==0:
            lsd = None
        else:
            lsd = max(int(round(abs(np.log10(abs(dvarmask).min()))+0.5))+6,6) #get decimal place of smallest value
            if abs(dvarmask).min()>=1.0:
                lsd=6
        spectra = ncd.createVariable("spectra","f8",("time","wvl"),
                                     zlib=True,least_significant_digit=lsd)
        spectra.set_auto_mask(False)
        spectra[:] = datavar[:]
        spectra.units = "erg cm-2 s-1 Hz-1"
        spectra.standard_name = "image_spectrum_avg"
        spectra.long_name = "image_spectrum_avg"
        
    else: #transits
        _log(logfile,"Writing a spectrum for each terminator column for each observation.....")
        datavar = dataset["transits"]
        dvarmask = datavar[np.isfinite(datavar)]
        dvarmask = dvarmask[dvarmask!=0]
        if len(dvarmask)==0:
            lsd = None
        else:
            lsd = max(int(round(abs(np.log10(abs(dvarmask).min()))+0.5))+6,6) #get decimal place of smallest value
            if abs(dvarmask).min()>=1.0:
                lsd=6
        transits = ncd.createVariable("transits","f8",("time","coord","wvl"),
                                      zlib=True,least_significant_digit=lsd)
        transits.set_auto_mask(False)
        transits[:] = datavar[:]
        transits.units = "km"
        transits.standard_name = "transit_spectra_map"
        transits.long_name = "transit_spectra_map"
        
        _log(logfile,"Recording the angular width of each column in degrees....")
        weights = ncd.createVariable("weights","f4",("time","coord"),
                                     zlib=True,least_significant_digit=6)
        weights.set_auto_mask(False)
        weights[:] = dataset["weights"][:]
        weights.units = "degrees"
        weights.standard_name = "transit_column_width"
        weights.long_name = "transit_column_width"
        
        _log(logfile,"Writing the terminator-averaged transit spectrum....")
        datavar = dataset["spectra"]
        dvarmask = datavar[np.isfinite(datavar)]
        dvarmask = dvarmask[dvarmask!=0]
        if len(dvarmask)==0:
            lsd = None
        else:
            lsd = max(int(round(abs(np.log10(abs(dvarmask).min()))+0.5))+6,6) #get decimal place of smallest value
            if abs(dvarmask).min()>=1.0:
                lsd=6
        spectra = ncd.createVariable("spectra","f8",("time","wvl"),
                                     zlib=True,least_significant_digit=lsd)
        spectra.set_auto_mask(False)
        spectra[:] = datavar[:]
        spectra.units = "km"
        spectra.standard_name = "transit_spectrum_avg"
        spectra.long_name = "transit_spectrum_avg"
    
    _log(logfile,"Write completed; handing off to parent routine.")
    ncd.sync()
    return ncd 
        
def _hdf5(filename,dataset,logfile=None,append=False):
    '''foo'''
    import h5py
    
    latitude  = dataset["lat"]
    longitude = dataset["lon"]
    time      = dataset["time"]
    wvls      = dataset["wvl"]
    
    _log(logfile,"Writing petitRADTRANS output to %s."%filename)
    _log(logfile,"Initializing file.....")
    
    mode = "w"
    if append:
        mode = "a"
    hdfile = h5py.File(filename,mode)
    
    if "lat" not in hdfile:
        if "images" in dataset:
           hdfile.create_dataset("lat",data=latitude[:].astype("float32"),compression='gzip',
                                 compression_opts=9,shuffle=True,fletcher32=True)
        else:
           hdfile.create_dataset("lat",data=latitude[:].astype("float32"),compression='gzip',
                                 maxshape=[None,dataset["lat"].shape[1]],compression_opts=9,
                                 shuffle=True,fletcher32=True)
        hdfile.attrs["lat"] = np.array(["latitude"       ,"degrees"])
    elif "transits" in hdfile:
        hdfile["lat"].resize(hdfile["lat"].shape[0]+dataset["lat"].shape[0],axis=0)
        hdfile["lat"][-dataset["lat"].shape[0]:] = dataset["lat"].astype("float32")
    if "lon" not in hdfile:
        if "images" in dataset:
           hdfile.create_dataset("lon",data=longitude[:].astype("float32"),compression='gzip',
                                 compression_opts=9,shuffle=True,fletcher32=True)
        else:
           hdfile.create_dataset("lon",data=longitude[:].astype("float32"),compression='gzip',
                                 maxshape=[None,dataset["lon"].shape[1]],compression_opts=9,
                                 shuffle=True,fletcher32=True)
        hdfile.attrs["lon"] = np.array(["longitude"      ,"degrees"])
    elif "transits" in hdfile:
        hdfile["lon"].resize(hdfile["lon"].shape[0]+dataset["lon"].shape[0],axis=0)
        hdfile["lon"][-dataset["lon"].shape[0]:] = dataset["lon"].astype("float32")
    if "wvl" not in hdfile:
        hdfile.create_dataset("wvl",data=wvls[:].astype("float32"),compression='gzip',
                            compression_opts=9,shuffle=True,fletcher32=True)
        hdfile.attrs["wvl"] = np.array(["wavelength"     ,"microns"])
    if "time" not in hdfile:
        hdfile.create_dataset("time",data=time[:].astype("float32"),compression='gzip',
                            maxshape=(None,),compression_opts=9,
                            shuffle=True,fletcher32=True)
        hdfile.attrs["time"]= np.array(["obsv_time_index","indices"])
    else:
        hdfile["time"].resize((hdfile["time"].shape[0]+len(time)),axis=0)
        hdfile["time"][-len(time):] = time[:]
    
    keyvars = list(dataset.keys())
    keyvars.remove("lat")
    keyvars.remove("lon")
    keyvars.remove("time")
    keyvars.remove("wvl")
    
    _log(logfile,"Packed fixed dimensions into %s\t...."%filename)
    
    for var in keyvars:
        if var not in hdfile:
            maxshape = [None,]
            for dim in dataset[var].shape[1:]:
                maxshape.append(dim)
            maxshape = tuple(maxshape)
            hdfile.create_dataset(var,data=dataset[var].astype("float64"),compression='gzip',
                                  maxshape=maxshape,compression_opts=9,shuffle=True,
                                  fletcher32=True)
            if "images" in keyvars and var=="spectra":
                hdfile.attrs[var] = np.array(["image_spectrum_avg","erg cm-2 s-1 Hz-1"])
            elif "transits" in keyvars and var=="spectra":
                hdfile.attrs[var] = np.array(["transit_spectrum_avg","km"])
            elif var=="images":
                hdfile.attrs[var] = np.array(["image_spectra_map","erg cm-2 s-1 Hz-1"])
            elif var=="colors":
                hdfile.attrs[var] = np.array(["true_color_image","RGB"])
            elif var=="transits":
                hdfile.attrs[var] = np.array(["transit_spectra_map","km"])
            elif var=="weights":
                hdfile.attrs[var] = np.array(["transit_column_width","degrees"])
        else:
            hdfile[var].resize((hdfile[var].shape[0]+dataset[var].shape[0]),axis=0)
            hdfile[var][-dataset[var].shape[0]:] = dataset[var].astype("float64")
        _log(logfile,"Packed %21s in %s\t...... %d timestamps"%(hdfile.attrs[var][0],
                                                                filename,dataset[var].shape[0])) 
    
    return hdfile

def save(filename,dataset,logfile=None,extracompression=False):
    '''Save petitRADTRANS ExoPlaSim output to a file.
        
    Output format is determined by the file extension in filename. Current supported formats are 
    NetCDF (\*.nc), HDF5 (\*.hdf5, \*.he5, \*.h5), numpy's ``np.savez_compressed`` format (\*.npz), and CSV format. If NumPy's 
    single-array .npy extension is used, .npz will be substituted--this is a compressed ZIP archive 
    containing .npy files. Additionally, the CSV output format can be used in compressed form either
    individually by using the .gz file extension, or collectively via tarballs (compressed or uncompressed).
    
    If a tarball format (e.g. \*.tar or \*.tar.gz) is used, output files will be packed into a tarball.
    gzip (.gz), bzip2 (.bz2), and lzma (.xz) compression types are supported. If a tarball format is 
    not used, then accepted file extensions are .csv, .txt, or .gz. All three will produce a directory
    named following the filename pattern, with one file per variable in the directory. If the .gz extension
    is used, NumPy will compress each output file using gzip compression. 
    
    CSV-type files will only contain 2D
    variable information, so the first N-1 dimensions will be flattened. The original variable shape is
    included in the file header (prepended with a # character) as the first items in a comma-separated
    list, with the first non-dimension item given as the '|||' placeholder. On reading variables from these
    files, they should be reshaped according to these dimensions. This is true even in tarballs (which 
    contain CSV files).
    
    A T21 model output with 10 vertical levels, 12 output times, all supported variables in grid 
    mode,and no standard deviation computation will have the following sizes for each format:
    
        +----------------+-----------+
        |Format          | Size      |
        +================+===========+
        |netCDF          | 12.8 MiB  |
        +----------------+-----------+
        |HDF5            | 17.2 MiB  |
        +----------------+-----------+
        |NumPy (default) | 19.3 MiB  |
        +----------------+-----------+
        |tar.xz          | 33.6 MiB  |
        +----------------+-----------+
        |tar.bz2         | 36.8 MiB  |
        +----------------+-----------+
        |gzipped         | 45.9 MiB  |
        +----------------+-----------+
        |uncompressed    | 160.2 MiB |
        +----------------+-----------+
            
    Using the NetCDF (.nc) format requires the netCDF4 python package.
    
    Using the HDF5 format (.h5, .hdf5, .he5) requires the h5py python package.
    
    Parameters
    ----------
    filename : str
        Path to the destination output file. The file extension determines the format. Currently,
        netCDF (\*.nc). numpy compressed (\*.npz), HDF5 (\*.hdf5, \*.he5, \*.h5), or CSV-type (\*.csv, \*.txt, \*.gz, \*.tar, \*.tar.gz,
        \*.tar.bz2, \*.tar.xz) are supported. If a format (such as npz) that requires
        that metadata be placed in a separate file is chosen, a second file with a '_metadata' suffix will be
        created.
    dataset : dict
        A dictionary containing the fields that should be written to output.
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
        
    Returns
    -------
    gcmt._Dataset object
        Open cross-format dataset object
    '''
    
    fileparts = filename.split('.')
    if fileparts[-1] == "nc":
        ncd = _netcdf(filename,dataset,logfile=logfile)
        ncd.close()
        return filename
    elif fileparts[-1] == "npz" or fileparts[-1] == "npy":
        meta,fn = _npsavez(filename,dataset,logfile=logfile)
        return fn
    elif (fileparts[-1] in ("csv","txt","gz","tar") or \
          (fileparts[-2]+"."+fileparts[-1]) in ("tar.gz","tar.bz2","tar.xz")):
        output = _csv(filename,dataset,logfile=logfile,extracompression=extracompression)
        return output
    elif fileparts[-1] in ("hdf5","h5","he5"):
        hdfile = _hdf5(filename,dataset,logfile=logfile)
        hdfile.close()
        return filename
    else:
        raise Exception("Unsupported output format detected. Supported formats are:\n\t\n\t%s"%("\n\t".join(pyburn.SUPPORTED)))
    
    
    
    
        
        
        
        
        
        
        
    
