import os
import argparse
from astropy.table import Table, vstack
from lssimpute import dirs, cat, impute
from LSS.tabulated_cosmo import TabulatedDESI
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser()
parser.add_argument('--tracer', '-t', default='LRG', help='Which LSS tracer catalogs to run on (LRG, ELG, QSO)')
parser.add_argument('--survey', '-s', default='y1mock', help='Survey to use (typically main or y1mock)')
parser.add_argument('--version', '-v', default=0, help='catalog version, for mocks this is mock number')
parser.add_argument('--impversion', '-i', default=None, help='override version for imputation, default None == same as version')
parser.add_argument('--nobackground', '-nb', action='store_true', help='Skip imputing "background" (close to random) galaxies.')
parser.add_argument('--overwrite', '-o', action='store_true', help='Set flag to allow overwriting of existing files')
parser.add_argument('--physical', '-p', action='store_true', help='Set flag to use physical units (S_perp, R) instead of Z and angular distance.')
parser.add_argument('--fit', '-f', action='store_true', help='Set flag to use model fit rather than KDE (physical units only).')
parser.add_argument('--radial_bins', '-rb', default=15, type=int, help='Number of radial bins to use')
parser.add_argument('--perp_bins', '-pb', default=18, type=int, help='Number of perpendicular bins to use')

# add dir management
# catdir (for base catalogs for reading, no writing)
# temp dir (intermediate files like logging, nn cats, etc)
# outdir (for imputed cats)

uargs = parser.parse_args()
catdir = dirs.get_catdir(uargs.survey, uargs.version)
if uargs.impversion is not None:
    impver = uargs.impversion
else:
    impver = uargs.version
stagedir = dirs.get_stagedir(uargs.survey, impver)
impute_dir = dirs.get_catdir('y1model', impver)

#read catalogs (using complete and full catalogs then splitting later)
clusn = Table.read(os.path.join(catdir, f'{uargs.tracer}_complete_N_clustering.dat.fits'))
cluss = Table.read(os.path.join(catdir, f'{uargs.tracer}_complete_S_clustering.dat.fits'))
full =  Table.read(os.path.join(catdir, f'{uargs.tracer}_full.dat.fits'))

#get missed and observed for each region
obsn, missn = cat.split_obs_missed(full, clusn, region='N')
obss, misss = cat.split_obs_missed(full, cluss, region='S')

#get nn cats, currently only extracting nn, not some set of them
mis_nncat_n = cat.extract_nn(missn, obsn)
mis_nncat_s = cat.extract_nn(misss, obss)
obs_nncat_n = cat.extract_nn(obsn)
obs_nncat_s = cat.extract_nn(obss)

# Add physical units
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

#merge both obs cause why not (maybe not idk for now)
obs_nncat = vstack([obs_nncat_n, obs_nncat_s])

#run imputation
impn = impute.ImputeModel(obs_nncat, mis_nncat_n)
impn_cat = impn.run(skip_background=uargs.nobackground, physical=uargs.physical, fit=uargs.fit, rbins=uargs.radial_bins, angbins=uargs.perp_bins)
figs = impn.figs
filename = f'{stagedir}/{uargs.tracer}_{uargs.survey}_{uargs.version}_N_model_bins_live.pdf'
with PdfPages(filename) as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close(fig)
imps = impute.ImputeModel(obs_nncat, mis_nncat_s)
imps_cat = imps.run(skip_background=uargs.nobackground, physical=uargs.physical, fit=uargs.fit, rbins=uargs.radial_bins, angbins=uargs.perp_bins)
figs = imps.figs
filename = f'{stagedir}/{uargs.tracer}_{uargs.survey}_{uargs.version}_S_model_bins_live.pdf'
with PdfPages(filename) as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close(fig)

imps.impute_details.write(os.path.join(stagedir, f'{uargs.tracer}_S_impute_details.fits'), overwrite=uargs.overwrite)
impn.impute_details.write(os.path.join(stagedir, f'{uargs.tracer}_N_impute_details.fits'), overwrite=uargs.overwrite)
imps_cat.write(os.path.join(impute_dir, f'{uargs.tracer}_S_clustering.dat.fits'), overwrite=uargs.overwrite)
impn_cat.write(os.path.join(impute_dir, f'{uargs.tracer}_N_clustering.dat.fits'), overwrite=uargs.overwrite)
