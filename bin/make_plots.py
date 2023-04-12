import os
import argparse
from astropy.table import Table, vstack
from lssimpute import plotting, dirs
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser()
parser.add_argument('--tracer', '-t', default='LRG', help='Which LSS tracer catalogs to run on (LRG, ELG, QSO)')
parser.add_argument('--survey', '-s', default='y1mock', help='Survey to use (typically main or y1mock)')
parser.add_argument('--version', '-v', default=0, help='catalog version, for mocks this is mock number')
parser.add_argument('--impversion', '-i', default=None, help='override version for imputation, default None == same as version')
parser.add_argument('--physical', '-p', action='store_true', help='Set flag to use physical units (S_perp, R) instead of Z and angular distance.')

uargs = parser.parse_args()
catdir = dirs.get_catdir(uargs.survey, uargs.version)
stagedir = dirs.get_stagedir(uargs.survey, uargs.version)
if uargs.impversion is not None:
    impver = uargs.impversion
else:
    impver = uargs.version
impute_dir = dirs.get_catdir('y1model', impver)

#read catalogs (using complete and full catalogs then splitting later)
clusn = Table.read(os.path.join(stagedir, f'{uargs.tracer}_N_clustering.dat.fits'))
cluss = Table.read(os.path.join(stagedir, f'{uargs.tracer}_S_clustering.dat.fits'))
missn = Table.read(os.path.join(stagedir, f'{uargs.tracer}_N_missing.dat.fits'))
misss = Table.read(os.path.join(stagedir, f'{uargs.tracer}_S_missing.dat.fits'))
full =  Table.read(os.path.join(catdir, f'{uargs.tracer}_full.dat.fits'))
impn = Table.read(os.path.join(impute_dir, f'{uargs.tracer}_N_clustering.dat.fits'))
imps = Table.read(os.path.join(impute_dir, f'{uargs.tracer}_S_clustering.dat.fits'))
detailsn = Table.read(os.path.join(stagedir, f'{uargs.tracer}_N_impute_details.fits'))
detailss = Table.read(os.path.join(stagedir, f'{uargs.tracer}_S_impute_details.fits'))

surv = f'{uargs.survey}/{uargs.version}'

plotter = plotting.validation_plots(imputecat=impn, cluscat=clusn, fullcat=full, miscat=missn, imputedetails=detailsn, survey=surv, tracer=uargs.tracer, region='N')
fig = plotter.imp_vs_true()
fig.savefig(f'{stagedir}/{uargs.tracer}_{uargs.survey}_{uargs.version}_N_impvtrue.pdf')
fig.close()
fig = plotter.n_z()
fig.savefig(f'{stagedir}/{uargs.tracer}_{uargs.survey}_{uargs.version}_N_nz.pdf')
plt.close(fig)
fig = plotter.fraction_bin()
fig.savefig(f'{stagedir}/{uargs.tracer}_{uargs.survey}_{uargs.version}_N_clusfrac.pdf')
plt.close(fig)
mode = 'physical' if uargs.physical else None
figs = plotter.imputation_fits(mode=mode)
filename = f'{stagedir}/{uargs.tracer}_{uargs.survey}_{uargs.version}_N_model_bins.pdf'
with PdfPages(filename) as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close(fig)

plotter = plotting.validation_plots(imputecat=imps, cluscat=cluss, fullcat=full, miscat=misss, imputedetails=detailss, survey=surv, tracer=uargs.tracer, region='S')
fig = plotter.imp_vs_true()
fig.savefig(f'{stagedir}/{uargs.tracer}_{uargs.survey}_{uargs.version}_S_impvtrue.pdf')
fig.close()
fig = plotter.n_z()
fig.savefig(f'{stagedir}/{uargs.tracer}_{uargs.survey}_{uargs.version}_S_nz.pdf')
plt.close(fig)
fig = plotter.fraction_bin()
fig.savefig(f'{stagedir}/{uargs.tracer}_{uargs.survey}_{uargs.version}_S_clusfrac.pdf')
plt.close(fig)
mode = 'physical' if uargs.physical else None
figs = plotter.imputation_fits(mode=mode)
filename = f'{stagedir}/{uargs.tracer}_{uargs.survey}_{uargs.version}_S_model_bins.pdf'
with PdfPages(filename) as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close(fig)