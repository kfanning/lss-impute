import os
import argparse
import numpy as np
from astropy.table import Table, vstack
from lssimpute import dirs, cat

parser = argparse.ArgumentParser()
parser.add_argument('--tracer', '-t', default='LRG', help='Which LSS tracer catalogs to run on (LRG, ELG, QSO)')
parser.add_argument('--survey', '-s', default='y1mock', help='Survey to use (typically main or y1mock)')
parser.add_argument('--version', '-v', default=0, help='catalog version, for mocks this is mock number')
parser.add_argument('--impversion', '-i', default=None, help='override version for imputation, default None == same as version')
parser.add_argument('--overwrite', '-o', action='store_true', help='Set flag to allow overwriting of existing files')
parser.add_argument('--stiched', action='store_true', help='Set flag to store randoms in stitched directory like stitch_impute.py')

# add dir management
# catdir (for base catalogs for reading, no writing)
# temp dir (intermediate files like logging, nn cats, etc)
# outdir (for imputed cats)

uargs = parser.parse_args()
catdir = dirs.get_catdir(uargs.survey, uargs.version)
stagedir = dirs.get_stagedir(uargs.survey, uargs.version)
if uargs.impversion is not None:
    impver = uargs.impversion
else:
    impver = uargs.version
impute_dir = dirs.get_catdir('y1model', impver)
if uargs.stiched:
    stitch_dir = dirs.get_catdir('y1impute', impver)
else:
    stiched_dir = impute_dir

#read catalogs (using impute and full to calc downsample amount)
impn = Table.read(os.path.join(impute_dir, f'{uargs.tracer}_N_clustering.dat.fits'))
imps = Table.read(os.path.join(impute_dir, f'{uargs.tracer}_S_clustering.dat.fits'))
full =  Table.read(os.path.join(catdir, f'{uargs.tracer}_full.dat.fits'))

# -1 Z in impute means NO imputation
imptgtids = np.append(impn[impn['Z']>0]['TARGETID'],imps[imps['Z']>0]['TARGETID'])
filled = full[(full['ZWARN']==0) | np.isin(full['TARGETID'], imptgtids)]
# get completeness/ntile
filled_ntiles = []
full_ntiles = []
for ntiles in range(1,8):
    filled_ntiles.append(len(filled[filled['NTILE']==ntiles]))
    full_ntiles.append(len(full[full['NTILE']==ntiles]))
filled_ntiles = np.array(filled_ntiles)
full_ntiles = np.array(full_ntiles)
comp_ntile = filled_ntiles/full_ntiles

# read random TODO add support for more than 1 random file
randcatn = Table.read(os.path.join(catdir, f'{uargs.tracer}_N_0_clustering.ran.fits'))
randcats = Table.read(os.path.join(catdir, f'{uargs.tracer}_S_0_clustering.ran.fits'))
new_randn = cat.downsample_randoms_by_ntile(randcatn, comp_ntile)
new_rands = cat.downsample_randoms_by_ntile(randcats, comp_ntile)
new_rands.write(os.path.join(stitch_dir, f'{uargs.tracer}_S_0_clustering.ran.fits'), format='fits', overwrite=uargs.overwrite)
new_randn.write(os.path.join(stitch_dir, f'{uargs.tracer}_N_0_clustering.ran.fits'), format='fits', overwrite=uargs.overwrite)
