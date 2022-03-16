import numpy as np
from scipy import interpolate

#From http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
widegamut = np.array([[ 1.4628067, -0.1840623, -0.2743606],
                      [-0.5217933,  1.4472381,  0.0677227],
                      [ 0.0349342, -0.0968930,  1.2884099]])




def _loadcie():
    wvl,xx,yy,zz = np.loadtxt("/".join(__file__.split("/")[:-1])+"/cmf.csv",unpack=True,delimiter=',')
    return wvl,xx,yy,zz
    
def interpolant(wvl):
    '''Construct interpolated color-matching functions for a series of wavelengths.
    
    Parameters
    ----------
    wvl : array-like
        Array of N wavelengths in nanometers
    
    Returns
    -------
    (array-like,array-like,array-like)
        (fx(lambda),fy(lambda),fz(lambda))
    '''
    w0,xx,yy,zz = _loadcie()
    fx = interpolate.interp1d(w0,xx)
    fy = interpolate.interp1d(w0,yy)
    fz = interpolate.interp1d(w0,zz)
    imin = np.where(wvl>np.amin(w0))[0][0]
    imax = np.where(wvl<np.amax(w0))[0][-1]
    wn = wvl[imin:imax+1]
    xn = fx(wn)
    yn = fy(wn)
    zn = fz(wn)
    return (xn,yn,zn)
    
def makexyz(wvl,spec,interpolant=None):
    '''Convert a spectrum to XYZ colour coordinates.
    
    The XYZ colour coordinate system is related to how the human eye's colour-receptive cells respond to
    light. The XYZ coordinates are computed by convolving the spectrum with three different empirically-derived
    response functions. These coordinates can then be transformed into RGB colour coordinates. Note that in this 
    system, (x,y,z) are brightness-normalized, while (X,Y,Z) are not. Additionally, z=1-x-y. Therefore, the
    three coordinates that are needed to produce an RGB tuple are not (x,y,z), but (x,y,Y). 
    
    Parameters
    ----------
    wvl : array-like
        Wavelengths in nanometers, of shape (N,)
    spec : array-like
        Spectrum, of shape (N,); units are arbitrary, but it should be given in flux, not flux density.
    interpolant : array-like, optional
        2D array of shape (3,N), corresponding to interpolated color-matching functions
        
    Returns
    -------
    numpy.ndarray, numpy.ndarray, numpy.ndarray
        x,y,Y--z can be inferred from x and y (z = 1-x-y), but Y preserves intensity information.
    '''
    if np.amin(wvl)<1.0e-3: #probably meters not nanometers
        wvl*=1.0e9
    w0,xx,yy,zz = _loadcie()
    imin = np.where(wvl>np.amin(w0))[0][0]
    imax = np.where(wvl<np.amax(w0))[0][-1]
    wn = wvl[imin:imax+1]
    specn = spec[imin:imax+1]
    if not interpolant:
        fx = interpolate.interp1d(w0,xx)
        fy = interpolate.interp1d(w0,yy)
        fz = interpolate.interp1d(w0,zz)
        
        xn = fx(wn)
        yn = fy(wn)
        zn = fz(wn)
    else:
        xn = interpolant[0]
        yn = interpolant[1]
        zn = interpolant[2]
    
    XI = np.trapz(xn[~np.isnan(specn)]*specn[~np.isnan(specn)],x=wn[~np.isnan(specn)])
    YI = np.trapz(yn[~np.isnan(specn)]*specn[~np.isnan(specn)],x=wn[~np.isnan(specn)])
    ZI = np.trapz(zn[~np.isnan(specn)]*specn[~np.isnan(specn)],x=wn[~np.isnan(specn)])
    xyzmin = np.amin((XI,YI,ZI))
    if xyzmin<0:
        XI -=xyzmin
        YI -=xyzmin
        ZI -=xyzmin
    if (XI+YI+ZI)>0:
        xnu = XI/(XI+YI+ZI)
        ynu = YI/(XI+YI+ZI)
        znu = 1.0-(xnu+ynu)
    else:
        xnu=0
        ynu=0
        znu=0
    
    return xnu,ynu,YI

def xyz2rgb(x,y,normalization):
    '''Convert (x,y) coordinates to RGB tuples, normalized to a given value.
    
    Note that z=1-x-y. This routine uses a wide gamut colourspace found at http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html.
    
    Parameters
    ----------
    x : array-like
        x colour-coordinate
    y : array-like
        y colour-coordinate.
    
    Returns
    -------
    array-like, array-like, array-like
        R,G,B colour values.
    '''
    
    z = 1-(x+y)
    
    r = np.sum(widegamut[0,:]*x)*normalization
    g = np.sum(widegamut[1,:]*y)*normalization
    b = np.sum(widegamut[2,:]*z)*normalization
    cmax = np.amax((r,g,b))
    
    #return r/cmax,g/cmax,b/cmax
    return r,g,b


def spec2rgb(wvl,spec,normalization=None):
    '''Convert a spectrum to (R,G,B) tuple, with optional normalization
    
    Parameters
    ----------
    wvl : array-like
        Array of length N containing wavelengths in nanometers
    spec : array-like
        Array of length N containing fluxes
    normalization : float, optional
        Maximum value. If not specified, defaults to 1.0.
        
    Returns
    -------
    (float,float,float)
        (R,G,B).
    '''
    x,y,I = makexyz(wvl,spec)
    if normalization:
        norm = normalization
    else:
        norm = 1.0
    
    r,g,b = xyz2rgb(x,y,norm)
    return r,g,b


def specs2rgb(wvl,specs):
    '''Convert a set of spectra into RGB colour values, so that relative intensity is preserved.
    
    Parameters
    ----------
    wvl : array-like
        Wavelengths in nanometers. Must have shape (N,)
    specs : array-like
        Spectra (fluxes), of shape (M,N) where M is the number of spectra.
        
    Returns
    -------
    array-like
        (M,3)-shape numpy array of R/G/B values
    '''
    interpol = interpolant(wvl)
    intensities = np.zeros((len(specs),3))
    for n in range(0,len(specs)):
        intensities[n,:] = makexyz(wvl,specs[n,:],interpolant=interpol)
    norms = intensities[:,2]/np.amax(intensities[:,2])
    colors = np.zeros((len(specs),3))
    for n in range(0,len(specs)):
        colors[n,:] = xyz2rgb(intensities[n,0],intensities[n,1],norms[n])
    
    return colors
    
    