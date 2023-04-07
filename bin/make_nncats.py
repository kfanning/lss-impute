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

# add dir management
# catdir (for base catalogs for reading, no writing)
# temp dir (intermediate files like logging, nn cats, etc)
# outdir (for imputed cats)

uargs = parser.parse_args()
catdir = dirs.get_catdir(uargs.survey, uargs.version)
stagedir = dirs.get_stagedir(uargs.survey, uargs.version)

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
obs_nncat_n = cat.extract_nn(missn, obsn)
obs_nncat_s = cat.extract_nn(misss, obss)

cosmo = TabulatedDESI()

mis_nncat_n['R'] = cosmo.comoving_radial_distance(mis_nncat_n['Z'])
mis_nncat_n['r_n0'] = cosmo.comoving_radial_distance(mis_nncat_n['z_n0'])
mis_nncat_n['sperp_n0'] = mis_nncat_n['r_n0']*mis_nncat_n['angdist_n0'] #small angle aprox + angdist in rad!
obs_nncat_n['r_n0'] = cosmo.comoving_radial_distance(obs_nncat_n['z_n0'])
obs_nncat_n['sperp_n0'] = obs_nncat_n['r_n0']*obs_nncat_n['angdist_n0'] #small angle aprox + angdist in rad!
obs_nncat_n['R'] = cosmo.comoving_radial_distance(obs_nncat_n['Z'])

mis_nncat_s['R'] = cosmo.comoving_radial_distance(mis_nncat_s['Z'])
mis_nncat_s['r_n0'] = cosmo.comoving_radial_distance(mis_nncat_s['z_n0'])
mis_nncat_s['sperp_n0'] = mis_nncat_s['r_n0']*mis_nncat_s['angdist_n0'] #small angle aprox + angdist in rad!
obs_nncat_s['r_n0'] = cosmo.comoving_radial_distance(obs_nncat_s['z_n0'])
obs_nncat_s['sperp_n0'] = obs_nncat_s['r_n0']*obs_nncat_s['angdist_n0'] #small angle aprox + angdist in rad!
obs_nncat_s['R'] = cosmo.comoving_radial_distance(obs_nncat_s['Z'])

mis_nncat_n.write(os.path.join(stagedir, f'{uargs.tracer}_N_missing.dat.fits'), overwrite=uargs.overwrite)
mis_nncat_s.write(os.path.join(stagedir, f'{uargs.tracer}_S_missing.dat.fits'), overwrite=uargs.overwrite)
obs_nncat_n.write(os.path.join(stagedir, f'{uargs.tracer}_N_clustering.dat.fits'), overwrite=uargs.overwrite)
obs_nncat_s.write(os.path.join(stagedir, f'{uargs.tracer}_S_clustering.dat.fits'), overwrite=uargs.overwrite)
