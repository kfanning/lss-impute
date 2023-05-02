import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, join
from .impute import ImputeModel

def uselato(self, fdir):
    import matplotlib as mpl
    for font in mpl.font_manager.findSystemFonts(fdir):
        mpl.font_manager.fontManager.addfont(font)
    mpl.rcParams['font.sans-serif'] = 'Lato'
    return

class validation_plots():

    def __init__(self, imputecat=None, fullcat=None, cluscat=None, miscat=None, imputedetails=None, survey=None, tracer=None, region=None):
        # stuff - merge NS catlogs beforehand
        # cluscat shoulw be the "staging" IE with  nn info clus cat"
        self.dpi = 200
        self.colors = {'complete': 'k', 'observed': 'c', 'impute_only': 'r', 'impute': 'm', 'missing': 'b'}    
        self.names = ['complete', 'observed', 'impute_only', 'impute', 'missing']
        self.cats = {'complete': fullcat, 'observed': cluscat, 'impute_only': imputecat, 'impute': None, 'missing': miscat}
        if cluscat is not None and imputecat is not None:
            self.cats['impute'] = vstack([cluscat, imputecat])
        # Just get rid of empty cats here
        self.names = [name for name in self.names if self.cats[name] is not None]
        self.survey = survey
        self.tracer = tracer
        self.region = region
        self.imputedetails = imputedetails
        if self.imputedetails is not None:
            self.fit = False
            for col in self.imputedetails.columns:
                if 'FIT' in col:
                    self.fit = True
                    break
        return

    def n_z(self):
        fig, ax = plt.subplots()
        fig.dpi = self.dpi
        nbins = 50
        zlo, zhi = self._find_zlims()
        for name in self.names:
            if name == 'complete':
                col = 'Z_not4clus'
            else:
                col = 'Z'
            if 'impute' in name:
                mask = np.array(self.cats[name]['Z']) > 0.0
            else:
                mask = np.ones(len(self.cats[name]), dtype=bool)
            bincounts, edges = np.histogram(self.cats[name][mask][col], range=(zlo,zhi), bins=nbins, density=True)
            ax.hist(self.cats[name][mask][col], bins=edges, color=self.colors[name], label=name, histtype='step', density=True)
        ax.legend()
        ax.set_ylabel('normalized n(z)')
        ax.set_xlabel('z')
        #ax.set_yscale('log')
        #ax.set_ylim([-100,600])
        ax.set_title(f'{self.survey} {self.tracer} {self.region}')
        return fig

    def imp_vs_true(self):
        if 'missing' in self.names and 'impute_only' in self.names:
            #get full catalog in impute catalog
            cat = join(self.cats['missing'], self.cats['impute_only'], keys='TARGETID', table_names=['comp','imp'])
            fig, ax = plt.subplots()
            fig.dpi = self.dpi
            ax.set_title(f'{self.survey} {self.tracer} {self.region}')
            ax.set_xlabel('True vs Impute Z Diff')
            ax.set_ylabel('density')
            catcut = cat[np.isin(cat['TARGETID'], self.cats['impute_only']['TARGETID'])]
            catcut = catcut[catcut['Z_imp'] > 0.0]
            zdiff = catcut['Z_comp'] - catcut['Z_imp']
            ax.hist(zdiff, color=self.colors['impute_only'], histtype='step', density=True, bins=50)
            return fig
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

    def imputation_fits(self, catver='observed', mode='physical'):
        '''
        '''
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
        rdiffs = cat[rname] - cat[f'{rname.lower()}_n0']
        rmins = list(self.imputedetails[f'MIN_{rname}'])
        rmaxs = list(self.imputedetails[f'MAX_{rname}'])
        pmins = list(self.imputedetails[f'MIN_{perpcol}'])
        pmaxs = list(self.imputedetails[f'MAX_{perpcol}'])
        binnum = list(self.imputedetails[f'BIN_NUM'])
        if self.fit:
            amp = list(self.imputedetails['FIT_AMPLITUDE'])
            sig = list(self.imputedetails['FIT_SIGMA'])
            slope = list(self.imputedetails['FIT_SLOPE'])
            intercept = list(self.imputedetails['FIT_INTERCEPT'])
        figs = []
        #for j in range(len(pmins)):
        for i in range(len(rmins)):
            mask = (cat[rcatcol] < rmaxs[i]) & (cat[rcatcol] > rmins[i]) & (cat[pcatcol] > pmins[i]) & (cat[pcatcol] < pmaxs[i])
            rdiffs = cat[mask][rname] - cat[mask][f'{rname.lower()}_n0']
            fig, axs = plt.subplots(1,2)#, sharey=True)
            fig.dpi=self.dpi
            fig.suptitle(f'bin: {binnum[i]} / {rmins[i]:.3f}{runit} < {rname} < {rmaxs[i]:.3f}{runit}, {pmins[i]:.3f}{perpunit} < {perpname} < {pmaxs[i]:.3f}{perpunit}')
            axs[0].set_ylabel('fraction of galaxies in bin')
            axs[0].set_title('"Background" Pairs')
            
            clusmask = (rdiffs < backg) & (rdiffs > -1*backg)
            clus = rdiffs[clusmask]
            back = rdiffs[~clusmask]
            cbbins, cbedges = np.histogram(back, bins=50, density=False)
            ccbins, ccedges = np.histogram(clus, bins=50, density=False)
            y1 = cbbins
            y2 = ccbins
            if self.fit:
                x = np.linspace(np.min(clus), np.max(clus), 50)
                params = (amp[i], sig[i], slope[i], intercept[i])
                y = ImputeModel.model(x, params)
                axs[1].plot(x,y, 'k--', label='fit')
            else:
                y1 /= len(rdiffs)
                y1 /= len(rdiffs)

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
            if name == 'complete':
                col = 'Z_not4clus'
            else:
                col = 'Z'
            maxes.append(np.nanmax(self.cats[name][col]))
            mins.append(np.nanmin(self.cats[name][col]))
        return np.min(mins), np.max(maxes)