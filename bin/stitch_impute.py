import os
import argparse
from astropy.table import Table, vstack
from lssimpute import dirs, cat
from LSS.tabulated_cosmo import TabulatedDESI

parser = argparse.ArgumentParser()
parser.add_argument('--tracer', '-t', default='LRG', help='Which LSS tracer catalogs to run on (LRG, ELG, QSO)')
parser.add_argument('--survey', '-s', default='y1mock', help='Survey to use (typically main or y1mock)')
parser.add_argument('--version', '-v', default=0, help='catalog version, for mocks this is mock number')
parser.add_argument('--overwrite', '-o', action='store_true', help='Set flag to allow overwriting of existing files')
parser.add_argument('--impversion', '-i', default=None, help='override version for imputation, default None == same as version')
parser.add_argument('--cut', '-c', action='store_true', help='Set flag to cut stiched together catalog to redshift range')

# add dir management
# catdir (for base catalogs for reading, no writing)
# temp dir (intermediate files like logging, nn cats, etc)
# outdir (for imputed cats)

zrange = {'LRG':[0.4,1.1], 'ELG':[0.8,1.6]}

uargs = parser.parse_args()
catdir = dirs.get_catdir(uargs.survey, uargs.version)
if uargs.impversion is not None:
    impver = uargs.impversion
else:
    impver = uargs.version
stagedir = dirs.get_stagedir(uargs.survey, impver)
impute_dir = dirs.get_catdir('y1model', impver)
stitch_dir = dirs.get_catdir('y1impute', impver)

#read catalogs (using complete and full catalogs then splitting later)
clusn = Table.read(os.path.join(catdir, f'{uargs.tracer}_N_clustering.dat.fits'))
cluss = Table.read(os.path.join(catdir, f'{uargs.tracer}_S_clustering.dat.fits'))
impn = Table.read(os.path.join(impute_dir, f'{uargs.tracer}_N_clustering.dat.fits'))
imps = Table.read(os.path.join(impute_dir, f'{uargs.tracer}_S_clustering.dat.fits'))

stitchedn = vstack([clusn, impn], join_type='inner')
if uargs.cut:
    stichedn = stichedn[(stitchedn['Z'] > zrange[uargs.tracer][0]) & (stitchedn['Z'] < zrange[uargs.tracer][1])]
stitchedn.write(os.path.join(stitch_dir, f'{uargs.tracer}_N_clustering.dat.fits'))
stitcheds = vstack([cluss, imps], join_type='inner')
if uargs.cut:
    sticheds = sticheds[(stitcheds['Z'] > zrange[uargs.tracer][0]) & (stitcheds['Z'] < zrange[uargs.tracer][1])]
stitcheds.write(os.path.join(stitch_dir, f'{uargs.tracer}_S_clustering.dat.fits'))
