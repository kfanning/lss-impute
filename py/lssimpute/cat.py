from astropy.table import Table, join
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
import numpy as np

def extract_nn(cat1, cat2=None, nmax=1):
    # when getting obs nn of missing, cat1 is missing cat2 is obs
    start = 1
    coord1 = SkyCoord(ra=cat1['RA']*u.degree, dec=cat1['DEC']*u.degree)
    if cat2 is None:
        coord2 = coord1
    else:
        coord2 = SkyCoord(ra=cat2['RA']*u.degree, dec=cat2['DEC']*u.degree)
    retcat = cat1.copy()
    shift = 1 if cat2 is None else 0 # if no cat 2 self compare and 1st nn is self
    if cat2 is None:
        cat2 = cat1 # this is by ref - so cleans up code below. cannot check for None now!
    store_tree = True if nmax > 1 else False
    for i in range(0+shift, nmax+shift):
        idx, d2d, d3d = match_coordinates_sky(coord1, coord2, nthneighbor=i+1, storekdtree=store_tree)
        #index cat by 0, nth=1 means first nearest
        retcat[f'ra_n{i-shift}'] = cat2['RA'][idx]
        retcat[f'dec_n{i-shift}'] = cat2['DEC'][idx]
        retcat[f'z_n{i-shift}'] = cat2['Z'][idx]
        retcat[f'tgtid_n{i-shift}'] = cat2['TARGETID'][idx]
        retcat[f'angdist_n{i-shift}'] = d2d.radian
    return retcat

def split_obs_missed(fullcat, compcat, region='N'):
    keep_comp = ['TARGETID', 'Z', 'WEIGHT', 'NZ', 'WEIGHT_FKP']
    fullreg = fullcat[fullcat['PHOTSYS']==region]
    fullreg = join(fullreg, compcat[keep_comp], keys='TARGETID', join_type='outer')
    missed = fullreg[fullreg['ZWARN'] != 0]
    observ = fullreg[fullreg['ZWARN'] == 0]
    return observ, missed

def downsample_mask(base_mask, fraction, seed=None):
    '''
    Randomly re-adds a fraction of the masked values

    Inputs:
        base_mask - np.array like with binary (True/False or 1/0) values
        fraction - float fraction of masked values (0) to unmask (change to 1)
    Out:
        new mask array (np.array) with fewer masked (0) values but same length
    '''
    masked_idx, _ = np.where(base_mask == 0)
    unmask_num = int(len(masked_idx) * fraction) #just get an int close to faction
    # setup rng shuffle to not require that source was preshuffled
    rng = np.random.default_rng(seed=seed)
    rng.shuffle(masked_idx)
    # unmask_num of len(masked_idx) will be unmasked
    unmask_idx = masked_idx[:unmask_num]
    base_mask[unmask_idx] = 1
    return base_mask

def downsample_randoms_by_ntile(randcat, fraction_by_ntile):
    '''
    Downsamples random catalog by given amount based on given faction
    to keep per NTILE value

    Inputs:
        randcat - astropy.Table like random galaxy catalog
        fraction_by_ntile - list like iterable with fraction of randoms to keep
                            per NTILE ordered from NTILE = 1 to NTILE = len()
    Out:
        Downsampled random catalog
    '''
    tot_mask = np.ones(len(randcat))
    for ntile, frac in enumerate(fraction_by_ntile, start=1):
        # mask out values at ntiles
        new_mask = np.array(randcat['NTILE']!=ntile)
        # downsample_mask removes frac of masked values
        downsampled = downsample_mask(new_mask, frac)
        # apply mask to intialized total mask
        tot_mask = tot_mask & downsampled
    return randcat[tot_mask] 

