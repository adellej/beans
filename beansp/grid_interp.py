from multiepoch_mcmc import grid_interpolator
import astropy.units as u
import astropy.constants as const
import numpy as np

c = const.c.to('cm s-1')
G = const.G.to('cm3 g-1 s-2')

def grid_interp(base, z, x_0, mdot, mass, radius, **kwargs):
    """
    Replacement for settle, based on Zac Johnston's Kepler model grid.
    Assumes a GridInterpolator object has been passed, keyword 'interpolator'
    
    Use the accompanying prior prior_grid with the ranges of the grid,
    which avoids evaluations outside the grid except for mdot
    
    :param base: base flux [MeV/nucleon] 
    :param z: accreted CNO metallicity
    :param x_0: accreted H-fraction
    :param mdot: accretion rate as a fraction of Eddington (as defined by the model code)
    :param mass: neutron star mass [M_sun]
    :param radius: neutron star radius [km]
    :param corr: correction function to modify the settle output (replaces the o
    
    :returns: 3-element array with recurrence time [hr], burst energy [1e39 erg], and alpha, all in the observer frame; OR None if the parameters fall outside the grid
    """

    if 'interpolator' not in kwargs:
        print ('** ERROR ** need to pass the interpolator to grid_interp')
        return None
    
    gi = kwargs['interpolator']
    if type(gi) != grid_interpolator.GridInterpolator:
        print ('** ERROR ** need a GridInterpolator object for the interpolation')
        return None
    
    # calculate redshift and g for the gridding
    
    R = radius*1e5*u.cm #cgs
    M = mass*const.M_sun.to('g') #cgs 
    redshift = np.power((1 - (2*G*M/(R*c**2))), -0.5).value
    gravity = (M*redshift*G/R**2 / (u.cm/u.s**2)).value / 1e14 #cgs

    try:
        # parameter ordering: ('mdot', 'qb', 'x', 'z', 'g')
        # print ([mdot, base, x_0, z, gravity])
        _result = gi.interpolate([mdot, base, x_0, z, gravity])[0]
        # result array ordering: ('rate', 'u_rate', 'energy', 'u_energy', 'peak', 'u_peak'))
    except ValueError:
        return None
    
    # convert here to observed parameters
    
    # and assemble result to return
    result = np.recarray(
        (1,), dtype=[("tdel", np.float64), ("E_b", np.float64), ("alpha", np.float64)]
    )
    # assign elements
    result.tdel = redshift*24./_result[0] # to hours, observer frame
    result.E_b = _result[2]/redshift/1e39 # to units of 1e39, observer frame
    result.alpha = 0.0 # don't currently calculate alpha
    result.mdot = mdot

    return result
