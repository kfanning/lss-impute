import os
import argparse
from astropy.table import Table, vstack
from lssimpute import dirs, cat, impute

parser = argparse.ArgumentParser()
parser.add_argument('--tracer', '-t', default='LRG', help='Which LSS tracer catalogs to run on (LRG, ELG, QSO)')
parser.add_argument('--survey', '-s', default='y1mock', help='Survey to use (typically main or y1mock)')
parser.add_argument('--version', '-v', default=0, help='catalog version, for mocks this is mock number')
parser.add_argument('--impversion', '-i', default=None, help='override version for imputation, default None == same as version')
parser.add_argument('--nobackground', '-nb', action='store_true', help='Skip imputing "background" (close to random) galaxies.')
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

#read catalogs (using complete and full catalogs then splitting later)
clusn = Table.read(os.path.join(catdir, f'{uargs.tracer}_complete_N_clustering.dat.fits'))
cluss = Table.read(os.path.join(catdir, f'{uargs.tracer}_complete_S_clustering.dat.fits'))
full =  Table.read(os.path.join(catdir, f'{uargs.tracer}_full.dat.fits'))

#get missed and observed for each region
obsn, missn = cat.split_obs_missed(full, clusn, region='N')
obss, misss = cat.split_obs_missed(full, cluss, region='S')

#merge both obs cause why not (maybe not idk for now)
obs = vstack([obsn, obss])

#get nn cats, currently only extracting nn, not some set of them
mis_nncat_n = cat.extract_nn(missn, obs)
mis_nncat_s = cat.extract_nn(misss, obs)
obs_nncat = cat.extract_nn(obs)

#run imputation
impn = impute.ImputeModel(obs_nncat, mis_nncat_n)
impn_cat = impn.run(skip_backbground=uargs.nobackground)
imps = impute.ImputeModel(obs_nncat, mis_nncat_s)
imps_cat = imps.run(skip_backbground=uargs.nobackground)

imps.impute_details.write(os.path.join(stagedir, f'{uargs.tracer}_S_impute_details.fits'))
impn.impute_details.write(os.path.join(stagedir, f'{uargs.tracer}_N_impute_details.fits'))
imps_cat.write(os.path.join(impute_dir, f'{uargs.tracer}_S_clustering.dat.fits'))
impn_cat.write(os.path.join(impute_dir, f'{uargs.tracer}_N_clustering.dat.fits'))
