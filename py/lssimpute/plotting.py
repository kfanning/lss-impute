import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

def uselato(self, fdir):
    import matplotlib as mpl
    for font in mpl.font_manager.findSystemFonts(fdir):
        mpl.font_manager.fontManager.addfont(font)
    mpl.rcParams['font.sans-serif'] = 'Lato'
    return

class validation_plots():

    def __init__(self, imputecat=None, fullcat=None, cluscat=None, imputedetails=None, survey=None, tracer=None, region=None):
        # stuff - merge NS catlogs beforehand
        # cluscat shoulw be the "staging" IE with  nn info clus cat"
        self.dpi = 200
        self.colors = {'full': 'k', 'clustering': 'c', 'impute': 'r'}    
        self.names = ['full', 'clustering', 'impute']
        self.cats = {'full': fullcat, 'clustering': cluscat, 'impute': imputecat}
        # Just get rid of empty cats here
        self.names = [name for name in self.names if self.cats[name] is not None]
        self.survey = survey
        self.tracer = tracer
        self.region = region
        self.imputedetails = imputedetails
        return

    def n_z(self):
        fig, ax = plt.subplots()
        fig.dpi = self.dpi
        nbins = 50
        zlo, zhi = self._find_zlims()
        for name in self.names:
            if name == 'full':
                col = 'Z_not4clus'
            else:
                col = 'Z'
            edges, bincounts = np.histogram(self.cats[name][col], range=(zlo,zhi), bins=nbins, density=True)
            ax.stairs(edges, bincounts, color=self.colors[name], label=name)
        ax.legend()
        ax.set_ylabel('normalized n(z)')
        ax.set_xlabel('z')
        #ax.set_yscale('log')
        #ax.set_ylim([-100,600])
        ax.set_title(f'{self.survey} {self.tracer} {self.region}')
        return fig

    def imp_vs_true(self):
        if 'full' in self.names and 'impute' in self.names:
            #get full catalog in impute catalog
            fig, ax = plt.subplots()
            fig.dpi = self.dpi
            ax.set_title(f'{self.survey} {self.tracer}')
        else:
            return

    def fraction_bin(self):
        if self.imputedetails is not None:
            fig, ax = plt.subplots()
            fig.dpi = self.dpi
            ax.plot(self.imputedetails['BIN_NUM'], self.imputedetails['CLUSTERED_FRAC'], 'b.')
            ax.set_xlabel('Bin Number')
            ax.set_ylabel('Fraction of "clustered" galaxies')
            ax.set_title(f'{self.survey} {self.tracer} {self.region}')
            return fig
        else:
            return

    def imputation_bins(self):
        return

    def imputation_fits(self, catver='clustering', mode='physical'):
        if mode == 'physical':
            rname = 'R'
            runit = 'Mpc/h'
            perpname = '$S_\perp$'
            perpunit = 'Mpc/h'
            perpcol = 'SPERPDIST'
            rcatcol = 'r_n0'
            pcatcol = 'sperp_n0'
            backg = 22.5 #fix this to read from impute class?
        else:
            rname = 'Z'
            runit = None
            perpname = 'Angular Dist'
            perpunit = 'rad'
            perpcol = 'ANGDIST'
            pcatcol = 'angdist_n0'
            rcatcol = 'z_n0'
            backg = 0.01
        cat = self.cats[catver]
        rdiffs = cat[rname] - cat[f'{rname.lower}_n0']
        rmins = list(self.imputedetails[f'MIN_{rname}'])
        rmaxs = list(self.imputedetails[f'MAX_{rname}'])
        pmins = list(self.imputedetails[f'MIN_{perpcol}'])
        pmaxs = list(self.imputedetails[f'MAX_{perpcol}'])
        binnum = list(self.imputedetails[f'MAX_{perpcol}'])
        figs = []
        for i in range(len(rmins)):
            for j in range(len(pmins)):
                mask = (cat[rcatcol] < rmaxs[i]) & (cat[rcatcol] > rmins[i]) & (cat[pcatcol] > pmins[j]) & (cat[pcatcol] < pmaxs[j])
                rdiffs = cat[mask][rname] - cat[mask][f'{rname.lower}_n0']
                fig, axs = plt.subplots(1,2)#, sharey=True)
                fig.dpi=self.dpi
                fig.suptitle(f'bin: {j+(i*(len(pmins)))} / {rmins[i]:.3f}{runit} < {rname} < {rmaxs[i]:.3f}{runit}, {pmins[j]:.3f}{perpunit} < {perpname} < {pmaxs[j]:.3f}{perpunit}')
                axs[0].set_ylabel('fraction of galaxies in bin')
                axs[0].set_title('"Background" Pairs')
                
                clusmask = (rdiffs < backg) & (rdiffs > -1*backg)
                clus = rdiffs[clusmask]
                back = rdiffs[~clusmask]
                cbbins, cbedges = np.histogram(rdiffs, bins=50)
                ccbins, ccedges = np.histogram(rdiffs, bins=50)
                y1 = cbbins/len(rdiffs)
                y2 = ccbins/len(rdiffs)

                axs[0].hist(cbedges[:-1], cbedges, weights=y1, color='b')
                #axs[0].hist(cbedges2[:-1], cbedges2, weights=np.split(y1, 2)[1], color='b')
                #axs[0].plot(x1, yfit1, 'k--')
                #axs[0].plot(bsample_x, bkde_y, 'r-')
                axs[0].set_xlabel(f'{rname} diff ({runit})')

                axs[1].set_title('"Clustered" Pairs')
                axs[1].hist(ccedges[:-1], ccedges, weights=y2, color='g')
                #axs[1].plot(x2, yfit2, 'k--')
                #if has_clustered:
                #    axs[1].plot(csample_x, ckde_y, 'r-')
                axs[1].set_xlabel(f'{rname} diff ({runit})')
                figs.append(fig)
        return figs

    def _find_zlims(self):
        maxes = []
        mins = []
        for name in self.names:
            if name == 'full':
                col = 'Z_not4clus'
            else:
                col = 'Z'
            maxes.append(np.nanmax(self.cats[name][col]))
            mins.append(np.nanmin(self.cats[name][col]))
        return np.min(mins), np.max(maxes)