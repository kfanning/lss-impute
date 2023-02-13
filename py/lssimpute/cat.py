from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
import numpy as np

def extract_nn(cat1, cat2=None, nmax=1):
    start = 1
    coord1 = SkyCoord(ra=cat1['RA']*u.degree, dec=cat1['DEC']*u.degree)
    if cat2 is None:
        coord2 = coord1
    else:
        coord2 = SkyCoord(ra=cat2['RA']*u.degree, dec=cat2['DEC']*u.degree)
    retcat = cat1.copy()
    shift = 1 if cat2 is None else 0 # if no cat 2 self compare and 1st nn is self
    for i in range(0+shift, nmax+shift):
        idx, d2d, d3d = match_coordinates_sky(coordm, coordc, nthneighbor=i+1, storekdtree=True)
        #index cat by 0, nth=1 means first nearest
        retcat[f'ra_n{i-shift}'] = cat2['RA'][idx]
        retcat[f'dec_n{i-shift}'] = cat2['DEC'][idx]
        retcat[f'z_n{i-shift}'] = cat2['Z'][idx]
        retcat[f'tgtid_n{i-shift}'] = cat2['TARGETID'][idx]
        retcat[f'angdist_n{i-shift}'] = d2d.radian
    return retcat

def split_obs_missed(fullcat, compcat, compcat2):
    keep_comp = ['TARGETID','Z']
    #split full then merge Z (for both cats even if obs has them)