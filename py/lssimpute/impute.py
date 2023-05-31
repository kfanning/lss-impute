import sys
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from LSS.tabulated_cosmo import TabulatedDESI
from scipy.stats import gaussian_kde
import scipy

def make_true():
    return

def make_nn_imp():
    return



class ImputeModel():

    def __init__(self, cluscat, misscat, region=None, version=None):
        self.cluscat = cluscat
        self.misscat = misscat
        self.cosmo = TabulatedDESI()
        self.misscat['r_n0'] = self.cosmo.comoving_radial_distance(self.misscat['z_n0'])
        self.misscat['sperp_n0'] = self.misscat['r_n0']*self.misscat['angdist_n0'] #small angle aprox + angdist in rad!
        self.cluscat['r_n0'] = self.cosmo.comoving_radial_distance(self.cluscat['z_n0'])
        self.cluscat['sperp_n0'] = self.cluscat['r_n0']*self.cluscat['angdist_n0'] #small angle aprox + angdist in rad!
        self.cluscat['R'] = self.cosmo.comoving_radial_distance(self.cluscat['Z'])
        self.did_angbin = False
        self.did_zbin = False
        self.did_sperpbin = False
        self.did_rbin = False
        self.backg_R = 22.5 #Mpc/h
        self.backg_Z = 0.01
        return

    def bin_angular(self, nbins=18):
        #Ang dist binning
        maxbin = np.max(self.misscat['angdist_n0'])#np.min((np.max(tab['angdist_n0']), np.max(clusN['angdist_n1'])))
        print(maxbin)
        #nbins = 18
        selectclus = self.cluscat[self.cluscat['angdist_n0'] < maxbin]
        self.ang_misbins, self.ang_edges = np.histogram(self.misscat['angdist_n0'], range=(0,maxbin), bins=nbins)
        self.ang_clusbins, clusedges = np.histogram(selectclus['angdist_n0'], range=(0,maxbin), bins=nbins)
        if nbins==18:
            # merge bins 6,7 and 8-14
            mergebin = 6
            maskb = np.ones(len(self.ang_misbins))
            #maskb[7] = 0
            maskb[mergebin:-1] = 0
            maske = np.ones(len(self.ang_edges))
            #maske[7] = 0
            maske[mergebin+1:-1] = 0
            #misbins[6] += misbins[7]
            #clusbins[6] += clusbins[7]
            self.ang_misbins[-1] += np.sum(self.ang_misbins[mergebin:])
            self.ang_clusbins[-1] += np.sum(self.ang_clusbins[mergebin:])
            self.ang_misbins = np.extract(maskb, self.ang_misbins)
            self.ang_edges = np.extract(maske, self.ang_edges)
            self.ang_clusbins = np.extract(maskb, self.ang_clusbins)
            self.did_angbin = True
        return

    def bin_sperp(self, nbins=18):
        #S_perp dist binning
        maxbin = np.max(self.misscat['sperp_n0'])
        #nbins = 18
        selectclus = self.cluscat[self.cluscat['sperp_n0'] < maxbin]
        self.sperp_misbins, self.sperp_edges = np.histogram(self.misscat['sperp_n0'], range=(0,maxbin), bins=nbins)
        self.sperp_clusbins, clusedges = np.histogram(selectclus['sperp_n0'], range=(0,maxbin), bins=nbins)
        if nbins==18:
            # merge bins 6,7 and 8-14
            mergebin = 11
            maskb = np.ones(len(self.sperp_misbins))
            #maskb[7] = 0
            maskb[mergebin:-1] = 0
            maske = np.ones(len(self.sperp_edges))
            #maske[7] = 0
            maske[mergebin+1:-1] = 0
            #misbins[6] += misbins[7]
            #clusbins[6] += clusbins[7]
            self.sperp_misbins[-1] += np.sum(self.sperp_misbins[mergebin:])
            self.sperp_clusbins[-1] += np.sum(self.sperp_clusbins[mergebin:])
            self.sperp_misbins = np.extract(maskb, self.sperp_misbins)
            self.sperp_edges = np.extract(maske, self.sperp_edges)
            self.sperp_clusbins = np.extract(maskb, self.sperp_clusbins)
        else:
            bins = np.linspace(0, 2. + 2./(nbins-1), nbins+1)
            bins[-1] = maxbin
            self.sperp_misbins, self.sperp_edges = np.histogram(self.misscat['sperp_n0'], bins=bins)
            self.sperp_clusbins, clusedges = np.histogram(selectclus['sperp_n0'], bins=bins)
        self.did_sperpbin = True
        return

    def plot_angbins(self, show=False):
        ntot = np.sum(self.ang_misbins) + np.sum(self.ang_clusbins)
        s_misbins = self.ang_misbins/ntot * 100
        s_clusbins = self.ang_clusbins/ntot * 100 #make percent
        #width = np.rad2deg(maxbin/nbins)*3600 #make arcsec
        width = np.rad2deg(self.ang_edges[1:]-self.ang_edges[:-1])*3600
        fig = plt.figure(dpi=200)
        plt.bar(np.rad2deg(self.ang_edges[:-1])*3600, s_clusbins, width=width, align='edge', color='b', label='Observed pairs')
        plt.bar(np.rad2deg(self.ang_edges[:-1])*3600, s_misbins, bottom=s_clusbins, width=width, align='edge', color='r', label='Pairs with 1 missed')
        plt.legend()
        plt.ylabel('Percent of nn pairs')
        plt.xlabel('nn pair distance (arcsec)')
        plt.title(f'Galaxies with pair distance < {np.rad2deg(max(self.ang_edges))*3600:.1f} arcsec \n No unobserved pairs above this')
        if show:
            plt.show()
        return fig

    def bin_z(self, nbins=15):
        #nn z binning
        # loooking at nn redshift
        maxbin = max([np.max(self.misscat['z_n0']), np.max(self.cluscat['z_n0'])]) #1.1
        minbin = min([np.min(self.misscat['z_n0']), np.min(self.cluscat['z_n0'])]) #0.4# might want to use range of all data?
        #nbins = 15#18 #20
        selectclus2 = self.cluscat[(self.cluscat['z_n0'] < maxbin) & (self.cluscat['z_n0'] > minbin)]
        self.z_misbins, self.z_edges = np.histogram(self.misscat['z_n0'], range=(minbin,maxbin), bins=nbins)
        self.z_clusbins, _ = np.histogram(selectclus2['z_n0'], range=(minbin,maxbin), bins=nbins)
        if nbins==15:
            # merge first + last bins
            num_merge = 3#4
            mask = np.ones(len(self.z_misbins))
            mask[1:num_merge+1] = 0
            mask[-1*num_merge-1:-1] = 0

            mask2 = np.ones(len(self.z_edges))
            mask2[1:num_merge+1] = 0
            mask2[-1*num_merge-1:-1] = 0

            self.z_misbins[0] = np.sum(self.z_misbins[0:num_merge+1])
            #misbins2[0] = misbins2[0] + misbins2[1] + misbins2[2] + misbins2[3]
            self.z_clusbins[0] = np.sum(self.z_clusbins[0:num_merge+1])
            #clusbins2[0] = clusbins2[0] + clusbins2[1] + clusbins2[2] + clusbins2[3]
            self.z_misbins[-1] = np.sum(self.z_misbins[-1*num_merge-1:-1])
            #misbins2[-1] = misbins2[-1] + misbins2[-2] + misbins2[-3] + misbins2[-4]
            self.z_clusbins[-1] = np.sum(self.z_clusbins[-1*num_merge-1:-1])
            #clusbins2[-1] = clusbins2[-1] + clusbins2[-2] + clusbins2[-3] + clusbins2[-4]
            #print(misbins2)
            self.z_misbins = np.extract(mask, self.z_misbins)
            self.z_clusbins = np.extract(mask, self.z_clusbins)
            self.z_edges = np.extract(mask2, self.z_edges)
        self.did_zbin = True
        return

    def bin_r(self, nbins=15):
        #nn r binning
        # loooking at nn r
        maxbin = max([np.max(self.misscat['r_n0']), np.max(self.cluscat['r_n0'])]) #1.1
        minbin = min([np.min(self.misscat['r_n0']), np.min(self.cluscat['r_n0'])]) #0.4# might want to use range of all data?
        #nbins = 15#18 #20
        selectclus2 = self.cluscat[(self.cluscat['r_n0'] < maxbin) & (self.cluscat['r_n0'] > minbin)]
        self.r_misbins, self.r_edges = np.histogram(self.misscat['r_n0'], range=(minbin,maxbin), bins=nbins)
        self.r_clusbins, _ = np.histogram(selectclus2['r_n0'], range=(minbin,maxbin), bins=nbins)
        if nbins==15:
            # merge first + last bins
            num_merge = 3#4
            mask = np.ones(len(self.r_misbins))
            mask[1:num_merge+1] = 0
            mask[-1*num_merge-1:-1] = 0

            mask2 = np.ones(len(self.r_edges))
            mask2[1:num_merge+1] = 0
            mask2[-1*num_merge-1:-1] = 0

            self.r_misbins[0] = np.sum(self.r_misbins[0:num_merge+1])
            #misbins2[0] = misbins2[0] + misbins2[1] + misbins2[2] + misbins2[3]
            self.r_clusbins[0] = np.sum(self.r_clusbins[0:num_merge+1])
            #clusbins2[0] = clusbins2[0] + clusbins2[1] + clusbins2[2] + clusbins2[3]
            self.r_misbins[-1] = np.sum(self.r_misbins[-1*num_merge-1:-1])
            #misbins2[-1] = misbins2[-1] + misbins2[-2] + misbins2[-3] + misbins2[-4]
            self.r_clusbins[-1] = np.sum(self.r_clusbins[-1*num_merge-1:-1])
            #clusbins2[-1] = clusbins2[-1] + clusbins2[-2] + clusbins2[-3] + clusbins2[-4]
            #print(misbins2)
            self.r_misbins = np.extract(mask, self.r_misbins)
            self.r_clusbins = np.extract(mask, self.r_clusbins)
            self.r_edges = np.extract(mask2, self.r_edges)
        self.did_rbin = True
        return

    def plot_zbins(self, show=False):
        ntot = np.sum(self.z_misbins) + np.sum(self.z_clusbins)
        s_misbins2 = self.z_misbins/ntot * 100
        s_clusbins2 = self.z_clusbins/ntot * 100 #make percent
        width = (max(self.z_edges)-max(self.z_edges))/len(self.z_edges)#np.rad2deg(maxbin/nbins)*3600 #make arcsec
        fig = plt.figure(dpi=200)
        plt.bar(self.z_edges[:-1], s_clusbins2, width=width, align='edge', color='b', label='Observed pairs')
        plt.bar(self.z_edges[:-1], s_misbins2, bottom=s_clusbins2, width=width, align='edge', color='r', label='Pairs with 1 missed')
        plt.legend()
        plt.ylabel('Percent of nn pairs')
        plt.xlabel('nn z')
        if show:
            plt.show()
        return fig 

    def impute(self, clusfrac_override=None, skip_background=False):
        #impute model
        mistab = self.misscat.copy()
        mistab['randnum'] = np.random.random_sample(len(mistab))
        mistab['Z'] = np.zeros(len(mistab), dtype=np.float64) - 1
        # data storage
        names = ['BIN_NUM', 'MIN_Z', 'MAX_Z', 'MIN_ANGDIST', 'MAX_ANGDIST', 'CLUSTERED_FRAC', 'N_OBS_CLUS', 'N_OBS_BACK', 'N_MIS_CLUS', 'N_MIS_BACK']
        binnum = []
        minzs = []
        maxzs = []
        minangs = []
        maxangs = []
        nominalfrac = []
        n_obsclus = []
        n_obsback = []
        n_misclus = []
        n_misback = []

        for j in range(len(self.ang_edges)-1): #misedges is nn_angdist
            for i in range(len(self.z_edges)-1): #misedges2 is nn_z
                maxz = self.z_edges[i+1]
                minz = self.z_edges[i]
                mindist = self.ang_edges[j]
                maxdist = self.ang_edges[j+1]
                # changed n1 to n0 since the nn extractor should shift for self matches now, be careful
                selclus = self.cluscat[(self.cluscat['z_n0'] > minz) & (self.cluscat['z_n0'] < maxz) & (self.cluscat['angdist_n0'] > mindist) & (self.cluscat['angdist_n0'] < maxdist)]# & (selectclus['Z'] > 0.4) & (selectclus['Z'] < 1.1)] #n1 is just cause the file has self as "n0" rather than a different object
                #selmis = mistab[(mistab['z_n0'] > minz) & (mistab['z_n0'] < maxz) & (mistab['angdist_n0'] > mindist) & (mistab['angdist_n0'] < maxdist)]# & (tab['RSDZ'] > 0.4) & (tab['RSDZ'] < 1.1)]

                selclus['zdiff'] = (selclus['Z']-selclus['z_n0'])
                #Sselmis['zdiff'] = (selmis['Z']-selmis['z_n0'])

                #Seperate cluster/background cutoff is just +/- 0.01 for now (by inspection :))
                backg = self.backg_Z# 0.01
                clusmask = (selclus['zdiff'] < backg) & (selclus['zdiff'] > -1*backg)
                clus_clus = selclus[clusmask]
                clus_back = selclus[~clusmask]

                #ccbins1, cbedges1 = np.histogram(clus_back['zdiff'], bins=25, range=(min(clus_back['zdiff']), -backg))
                #ccbins2, cbedges2 =  np.histogram(clus_back['zdiff'], bins=25, range=(backg, max(clus_back['zdiff'])))
                #cbbins = np.concatenate((ccbins1, ccbins2))
                #cbedges = np.concatenate((cbedges1, cbedges2))
                cbbins, cbedges = np.histogram(clus_back['zdiff'], bins=50)#, range=(-0.01,0.01))
                #cbbins1, cbedges1 = np.histogram(clus_back[clusback['zdiff'] < -1*backg]['zdiff'], bins=50)#, range=(-0.01,0.01))
                #cbbins2, cbedges2 = np.histogram(clus_back[clusback['zdiff'] > backg]['zdiff'], bins=50)#, range=(-0.01,0.01))
                ccbins, ccedges = np.histogram(clus_clus['zdiff'], bins=50)#, range=(-0.01,0.01))

                y1 = cbbins/(len(clus_back)*(cbedges[1]-cbedges[0]))
                x1 = ((cbedges[1:] + cbedges[:-1])/2)#np.concatenate()(cbedges2[1:] + cbedges2[:-1])/2))
                y2 = ccbins/(len(clus_clus)*(ccedges[1]-ccedges[0]))
                x2 = (ccedges[1:] + ccedges[:-1])/2

                bkde = gaussian_kde(clus_back['zdiff'])
                has_clustered = (len(clus_clus) > 1)
                if has_clustered:
                    ckde = gaussian_kde(clus_clus['zdiff'])#, h=0.01)
                    csample_x = np.linspace(-backg,backg,100)
                    ckde_y = ckde.evaluate(csample_x)
                else:
                    print('no "clustered" galaxies!')
                bsample_x = np.linspace(-1, 1, 100)
                bkde_y = bkde.evaluate(bsample_x)

                # Draw z's ###CHECK THIS####
                if clusfrac_override is not None:
                    clus_frac = clusfrac_override
                else:  
                    clus_frac = len(clus_clus)/len(selclus)

                #maxz = max(selectclus['Z'])
                #for row in fittab:
                mask = (mistab['z_n0'] < maxz) & (mistab['z_n0'] > minz) & (mistab['angdist_n0'] > mindist) & (mistab['angdist_n0'] < maxdist) & ((mistab['Z'] < 0))# | (mistab['Z'] > maxz)) #ensure positive
                while np.count_nonzero(mask) > 0:
                    mask = (mistab['z_n0'] < maxz) & (mistab['z_n0'] > minz) & (mistab['angdist_n0'] > mindist) & (mistab['angdist_n0'] < maxdist) & ((mistab['Z'] < 0))# | (mistab['Z'] > maxz)) #ensure positive
                    back = (mistab['randnum'] > clus_frac)
                    if has_clustered:
                        clus = (mistab['randnum'] < clus_frac)
                        mistab['Z'][mask & clus] = mistab[mask & clus]['z_n0'] + ckde.resample(np.count_nonzero(mask & clus))[0]
                    if skip_background:
                        mask = mask & (~back)
                    else:
                        mistab['Z'][mask & back] = mistab[mask & back]['z_n0'] + bkde.resample(np.count_nonzero(mask & back))[0]

                select_miss = mistab[(mistab['z_n0'] < maxz) & (mistab['z_n0'] > minz) & (mistab['angdist_n0'] > mindist) & (mistab['angdist_n0'] < maxdist)]
                zdiff_new = select_miss['Z'] - select_miss['z_n0']
                miss_clus_mask = (zdiff_new < backg) & (zdiff_new > -1*backg)
                # data collection
                binnum.append(i+(j*(len(self.z_edges)-1))) 
                minzs.append(minz)
                maxzs.append(maxz)
                minangs.append(mindist)
                maxangs.append(maxdist)
                nominalfrac.append(clus_frac)
                n_obsclus.append(len(clus_clus))
                n_obsback.append(len(clus_back))
                n_misclus.append(len(select_miss[miss_clus_mask]))
                n_misback.append(len(select_miss[~miss_clus_mask]))

                if False: #((i % 5) == 0) | (i > 0): 
                    fig, axs = plt.subplots(1,2)#, sharey=True)
                    fig.dpi=200
                    fig.suptitle(f'{minz:.3f} < Z < {maxz:.3f}, {np.rad2deg(mindist)*3600:.3f} arcsec < Angular Seperation < {np.rad2deg(maxdist)*3600:.3f} arcsec')

                    axs[0].set_ylabel('probability density of galaxies')
                    axs[0].set_title('"Background" Pairs')
                    axs[0].hist(cbedges[:-1], cbedges, weights=y1, color='b')
                    #axs[0].hist(cbedges2[:-1], cbedges2, weights=np.split(y1, 2)[1], color='b')
                    #axs[0].plot(x1, yfit1, 'k--')
                    axs[0].plot(bsample_x, bkde_y, 'r-')
                    axs[0].set_xlabel('Z diff')

                    axs[1].set_title('"Clustered" Pairs')
                    axs[1].hist(ccedges[:-1], ccedges, weights=y2, color='g')
                    #axs[1].plot(x2, yfit2, 'k--')
                    if has_clustered:
                        axs[1].plot(csample_x, ckde_y, 'r-')
                    axs[1].set_xlabel('Z diff')
        
        self.impute_details = Table([binnum, minzs, maxzs, minangs, maxangs, nominalfrac, n_obsclus, n_obsback, n_misclus, n_misback], names=names)
        #mistab.write('impute_model_20230224_5S.fits', format='fits', overwrite=True)
        #fittab.write('model_fit_20220609.fits', format='fits', overwrite=True)
        return mistab

    def impute_physical(self, clusfrac_override=None, skip_background=False):
        #impute model
        mistab = self.misscat.copy()
        mistab['randnum'] = np.random.random_sample(len(mistab))
        mistab['R'] = np.zeros(len(mistab), dtype=np.float64) - 1
        mistab['Z'] = np.zeros(len(mistab), dtype=np.float64) - 1
        # data storage
        names = ['BIN_NUM', 'MIN_R', 'MAX_R', 'MIN_SPERPDIST', 'MAX_SPERPDIST', 'CLUSTERED_FRAC', 'N_OBS_CLUS', 'N_OBS_BACK', 'N_MIS_CLUS', 'N_MIS_BACK']
        binnum = []
        minrs = []
        maxrs = []
        minsperps = []
        maxsperps = []
        nominalfrac = []
        n_obsclus = []
        n_obsback = []
        n_misclus = []
        n_misback = []

        for j in range(len(self.sperp_edges)-1): #misedges is nn_angdist
            for i in range(len(self.r_edges)-1): #misedges2 is nn_z
                maxr = self.r_edges[i+1]
                minr = self.r_edges[i]
                minsperp = self.sperp_edges[j]
                maxsperp = self.sperp_edges[j+1]
                # changed n1 to n0 since the nn extractor should shift for self matches now, be careful
                selclus = self.cluscat[(self.cluscat['r_n0'] > minr) & (self.cluscat['r_n0'] < maxr) & (self.cluscat['sperp_n0'] > minsperp) & (self.cluscat['sperp_n0'] < maxsperp)]# & (selectclus['Z'] > 0.4) & (selectclus['Z'] < 1.1)] #n1 is just cause the file has self as "n0" rather than a different object
                selclus['rdiff'] = (selclus['R']-selclus['r_n0'])

                #Seperate cluster/background cutoff is just +/- 0.01 for now (by inspection :))
                #backg = 0.01
                # 0.01 zdiff at low Z is ~ 22.5Mpc/h, ~12.5Mpc/h at high z where it maybe doesn't work well
                backg = self.backg_R#22.5 #Mpc/h
                clusmask = (selclus['rdiff'] < backg) & (selclus['rdiff'] > -1*backg)
                clus_clus = selclus[clusmask]
                clus_back = selclus[~clusmask]

                #ccbins1, cbedges1 = np.histogram(clus_back['zdiff'], bins=25, range=(min(clus_back['zdiff']), -backg))
                #ccbins2, cbedges2 =  np.histogram(clus_back['zdiff'], bins=25, range=(backg, max(clus_back['zdiff'])))
                #cbbins = np.concatenate((ccbins1, ccbins2))
                #cbedges = np.concatenate((cbedges1, cbedges2))
                cbbins, cbedges = np.histogram(clus_back['rdiff'], bins=50)#, range=(-0.01,0.01))
                #cbbins1, cbedges1 = np.histogram(clus_back[clusback['zdiff'] < -1*backg]['zdiff'], bins=50)#, range=(-0.01,0.01))
                #cbbins2, cbedges2 = np.histogram(clus_back[clusback['zdiff'] > backg]['zdiff'], bins=50)#, range=(-0.01,0.01))
                ccbins, ccedges = np.histogram(clus_clus['rdiff'], bins=50)#, range=(-0.01,0.01))

                y1 = cbbins/(len(clus_back)*(cbedges[1]-cbedges[0]))
                x1 = ((cbedges[1:] + cbedges[:-1])/2)#np.concatenate()(cbedges2[1:] + cbedges2[:-1])/2))
                y2 = ccbins/(len(clus_clus)*(ccedges[1]-ccedges[0]))
                x2 = (ccedges[1:] + ccedges[:-1])/2

                bkde = gaussian_kde(clus_back['rdiff'])
                has_clustered = (len(clus_clus) > 1)
                if has_clustered:
                    ckde = gaussian_kde(clus_clus['rdiff'])#, h=0.01)
                    csample_x = np.linspace(-backg,backg,100)
                    ckde_y = ckde.evaluate(csample_x)
                else:
                    print('no "clustered" galaxies!')
                bsample_x = np.linspace(-1, 1, 100)
                bkde_y = bkde.evaluate(bsample_x)

                # Draw z's ###CHECK THIS####
                if clusfrac_override is not None:
                    clus_frac = clusfrac_override
                else:  
                    clus_frac = len(clus_clus)/len(selclus)
                #maxz = max(selectclus['Z'])
                #for row in fittab:
                mask = (mistab['r_n0'] < maxr) & (mistab['r_n0'] > minr) & (mistab['sperp_n0'] > minsperp) & (mistab['sperp_n0'] < maxsperp) & ((mistab['R'] < 0))# | (mistab['Z'] > maxz)) #ensure positive
                while np.count_nonzero(mask) > 0:
                    mask = (mistab['r_n0'] < maxr) & (mistab['r_n0'] > minr) & (mistab['sperp_n0'] > minsperp) & (mistab['sperp_n0'] < maxsperp) & ((mistab['R'] < 0))# | (mistab['Z'] > maxz)) #ensure positive
                    back = (mistab['randnum'] > clus_frac)
                    if has_clustered:
                        clus = (mistab['randnum'] < clus_frac)
                        mistab['R'][mask & clus] = mistab[mask & clus]['r_n0'] + ckde.resample(np.count_nonzero(mask & clus))[0]
                    if skip_background:
                        mask = mask & (~back)
                        if not has_clustered:
                            break
                    else:
                        mistab['R'][mask & back] = mistab[mask & back]['r_n0'] + bkde.resample(np.count_nonzero(mask & back))[0]

                select_miss = mistab[(mistab['r_n0'] < maxr) & (mistab['r_n0'] > minr) & (mistab['sperp_n0'] > minsperp) & (mistab['sperp_n0'] < maxsperp)]
                rdiff_new = select_miss['R'] - select_miss['r_n0']
                miss_clus_mask = (rdiff_new < backg) & (rdiff_new > -1*backg)
                mistab['Z'] = self._inverse_comoving_radial_dist(mistab['R'])
                # data collection
                binnum.append(i+(j*(len(self.r_edges)-1))) 
                minrs.append(minr)
                maxrs.append(maxr)
                minsperps.append(minsperp)
                maxsperps.append(maxsperp)
                nominalfrac.append(clus_frac)
                n_obsclus.append(len(clus_clus))
                n_obsback.append(len(clus_back))
                n_misclus.append(len(select_miss[miss_clus_mask]))
                n_misback.append(len(select_miss[~miss_clus_mask]))

                if False: #((i % 5) == 0) | (i > 0): 
                    fig, axs = plt.subplots(1,2)#, sharey=True)
                    fig.dpi=200
                    fig.suptitle(f'{minr:.3f} < R < {maxr:.3f}, {minsperp} < S_perp < {maxsperp:.3f}')

                    axs[0].set_ylabel('probability density of galaxies')
                    axs[0].set_title('"Background" Pairs')
                    axs[0].hist(cbedges[:-1], cbedges, weights=y1, color='b')
                    #axs[0].hist(cbedges2[:-1], cbedges2, weights=np.split(y1, 2)[1], color='b')
                    #axs[0].plot(x1, yfit1, 'k--')
                    axs[0].plot(bsample_x, bkde_y, 'r-')
                    axs[0].set_xlabel('Z diff')

                    axs[1].set_title('"Clustered" Pairs')
                    axs[1].hist(ccedges[:-1], ccedges, weights=y2, color='g')
                    #axs[1].plot(x2, yfit2, 'k--')
                    if has_clustered:
                        axs[1].plot(csample_x, ckde_y, 'r-')
                    axs[1].set_xlabel('Z diff')
        
        self.impute_details = Table([binnum, minrs, maxrs, minsperps, maxsperps, nominalfrac, n_obsclus, n_obsback, n_misclus, n_misback], names=names)
        #mistab.write('impute_model_20230224_5S.fits', format='fits', overwrite=True)
        #fittab.write('model_fit_20220609.fits', format='fits', overwrite=True)
        return mistab

    def impute_physical_fit(self, clusfrac_override=None, skip_background=True, fit_type='gauss'):
        fit_type = fit_type.lower()
        assert fit_type in ['gauss', 'lorentz', 'quad', 'dynamic'], f'Invalid fit_type {fit_type}'

        #impute model
        mistab = self.misscat.copy()
        mistab['randnum'] = np.random.random_sample(len(mistab))
        mistab['R'] = np.zeros(len(mistab), dtype=np.float64) - 1
        mistab['Z'] = np.zeros(len(mistab), dtype=np.float64) - 1
        # data storage
        names = ['BIN_NUM', 'MIN_R', 'MAX_R', 'MIN_SPERPDIST', 'MAX_SPERPDIST', 'CLUSTERED_FRAC', 'CALC_FRAC', 'N_OBS_CLUS', 'N_OBS_BACK', 'N_MIS_CLUS', 'N_MIS_BACK', 'FIT_AMPLITUDE', 'FIT_WIDTH', 'FIT_SLOPE', 'FIT_INTERCEPT', 'FIT_QUAD', 'FIT_TYPE','ERROR']
        binnum = []
        minrs = []
        maxrs = []
        minsperps = []
        maxsperps = []
        nominalfrac = []
        fitfrac = []
        n_obsclus = []
        n_obsback = []
        n_misclus = []
        n_misback = []
        amp = []
        sig = []
        slope = []
        intercept = []
        quad = []
        fitt = []
        error_fit = []
        self.figs = []

        for j in range(len(self.sperp_edges)-1): #misedges is nn_angdist
            for i in range(len(self.r_edges)-1): #misedges2 is nn_z
                maxr = self.r_edges[i+1]
                minr = self.r_edges[i]
                minsperp = self.sperp_edges[j]
                maxsperp = self.sperp_edges[j+1]
                # changed n1 to n0 since the nn extractor should shift for self matches now, be careful
                selclus = self.cluscat[(self.cluscat['r_n0'] > minr) & (self.cluscat['r_n0'] < maxr) & (self.cluscat['sperp_n0'] > minsperp) & (self.cluscat['sperp_n0'] < maxsperp)]# & (selectclus['Z'] > 0.4) & (selectclus['Z'] < 1.1)] #n1 is just cause the file has self as "n0" rather than a different object
                selclus['rdiff'] = (selclus['R']-selclus['r_n0'])
                #Seperate cluster/background cutoff is just +/- 0.01 for now (by inspection :))
                #backg = 0.01
                # 0.01 zdiff at low Z is ~ 22.5Mpc/h, ~12.5Mpc/h at high z where it maybe doesn't work well
                backg = self.backg_R#22.5 #Mpc/h
                clusmask = (selclus['rdiff'] < backg) & (selclus['rdiff'] > -1*backg)
                clus_clus = selclus[clusmask]
                clus_back = selclus[~clusmask]

                cbbins, cbedges = np.histogram(clus_back['rdiff'], bins=50)#, range=(-0.01,0.01))
                ccbins, ccedges = np.histogram(clus_clus['rdiff'], bins=50)#, range=(-0.01,0.01))

                if fit_type=='dynamic':
                    ft = ['gauss','lorentz']#,'quad']
                else:
                    ft = [fit_type]
                errs = []
                ress = []
                for f in ft:
                    res, end_err = self._fit(clus_clus['rdiff'], fit_type=f)
                    errs.append(end_err*self.dof(f))
                    ress.append(res)
                    print(f'Fit: {f}, end error/dof: {end_err:3f}')
                    print(f'{i+(j*(len(self.r_edges)-1))} optimize status: {res.success}, {res.x}, {res.message}')
                choice_idx = np.argmin(errs)
                fit_choice = ft[choice_idx]
                res = ress[choice_idx]
                error_fit.append(errs[choice_idx]/self.dof(fit_choice))
                y1 = cbbins/(len(clus_back)*(cbedges[1]-cbedges[0]))
                x1 = ((cbedges[1:] + cbedges[:-1])/2)#np.concatenate()(cbedges2[1:] + cbedges2[:-1])/2))
                y2 = ccbins/(len(clus_clus)*(ccedges[1]-ccedges[0]))
                x2 = (ccedges[1:] + ccedges[:-1])/2

                has_clustered = (len(clus_clus) > 1)

                # Draw z's ###CHECK THIS####
                if clusfrac_override is not None:
                    clus_frac = clusfrac_override
                elif fit_choice in ['gauss','quad']:  
                    clus_frac = res.x[0]*res.x[1]*np.sqrt(2*np.pi)/len(selclus)
                elif fit_choice == 'lorentz':
                    clus_frac = res.x[0]*res.x[1]*np.pi/len(selclus)

                mask = (mistab['r_n0'] < maxr) & (mistab['r_n0'] > minr) & (mistab['sperp_n0'] > minsperp) & (mistab['sperp_n0'] < maxsperp) & ((mistab['R'] < 0))# | (mistab['Z'] > maxz)) #ensure positive
                while np.count_nonzero(mask) > 0:
                    mask = (mistab['r_n0'] < maxr) & (mistab['r_n0'] > minr) & (mistab['sperp_n0'] > minsperp) & (mistab['sperp_n0'] < maxsperp) & ((mistab['R'] < 0))# | (mistab['Z'] > maxz)) #ensure positive
                    back = (mistab['randnum'] > clus_frac)
                    if has_clustered:
                        clus = (mistab['randnum'] < clus_frac)
                        if fit_type in ['gauss','quad']:
                            sample = np.random.normal(0.0, res.x[1], np.count_nonzero(mask & clus))
                        elif fit_type == 'lorentz':
                            sample = res.x[1]*np.random.standard_cauchy(np.count_nonzero(mask & clus)) #lorentz works this way?
                        mistab['R'][mask & clus] = mistab[mask & clus]['r_n0'] + np.random.normal(0.0, res.x[1], np.count_nonzero(mask & clus))
                    if skip_background:
                        mask = mask & (~back)
                        if not has_clustered:
                            break
                
                select_miss = mistab[(mistab['r_n0'] < maxr) & (mistab['r_n0'] > minr) & (mistab['sperp_n0'] > minsperp) & (mistab['sperp_n0'] < maxsperp)]
                rdiff_new = select_miss['R'] - select_miss['r_n0']
                miss_clus_mask = (rdiff_new < backg) & (rdiff_new > -1*backg)
                mistab['Z'] = self._inverse_comoving_radial_dist(mistab['R'])

                #Plot
                if fit_choice == 'gauss':
                    model = self.model
                elif fit_choice == 'lorentz':
                    model = self.model_lorentz
                elif fit_choice == 'quad':
                    model = self.model_quad
                rname = 'R'
                runit = 'Mpc/h'
                perpname = '$S_\perp$'
                perpunit = 'Mpc/h'
                perpcol = 'SPERPDIST'
                rcatcol = 'r_n0'
                pcatcol = 'sperp_n0'
                fig, axs = plt.subplots(1,2)#, sharey=True)
                fig.dpi=200
                fig.suptitle(f'bin: {i+(j*(len(self.r_edges)-1))} / {minr:.3f}{runit} < {rname} < {maxr:.3f}{runit}, {minsperp:.3f}{perpunit} < {perpname} < {maxsperp:.3f}{perpunit}')
                axs[0].set_ylabel('fraction of galaxies in bin')
                axs[0].set_title('"Background" Pairs')
                x = np.linspace(np.min(clus_clus['rdiff']), np.max(clus_clus['rdiff']), 50)
                y = model(x, res.x)
                axs[1].plot(x,y, 'k-', label='fit')
                axs[0].hist(cbedges[:-1], cbedges, weights=cbbins, color='b')
                axs[0].set_xlabel(f'{rname} diff ({runit})')

                axs[1].set_title('"Clustered" Pairs')
                axs[1].hist(ccedges[:-1], ccedges, weights=ccbins, color='g')
                axs[1].set_xlabel(f'{rname} diff ({runit})')
                self.figs.append(fig)

                # data collection
                binnum.append(i+(j*(len(self.r_edges)-1))) 
                minrs.append(minr)
                maxrs.append(maxr)
                minsperps.append(minsperp)
                maxsperps.append(maxsperp)
                nominalfrac.append(len(clus_clus)/len(selclus))
                fitfrac.append(clus_frac)
                n_obsclus.append(len(clus_clus))
                n_obsback.append(len(clus_back))
                n_misclus.append(len(select_miss[miss_clus_mask]))
                n_misback.append(len(select_miss[~miss_clus_mask]))
                amp.append(res.x[0])
                sig.append(res.x[1])
                slope.append(res.x[2])
                intercept.append(res.x[3])
                if fit_choice == 'quad':
                    quad.append(res.x[4])
                else:
                    quad.append(0)
                fitt.append(fit_choice)
        self.impute_details = Table([binnum, minrs, maxrs, minsperps, maxsperps, nominalfrac, fitfrac, n_obsclus, n_obsback, n_misclus, n_misback, amp, sig, slope, intercept, quad, fitt, error_fit], names=names)
        #mistab.write('impute_model_20230224_5S.fits', format='fits', overwrite=True)
        #fittab.write('model_fit_20220609.fits', format='fits', overwrite=True)
        return mistab

    def _fit(self, data, fit_type='gauss'):
        #width = bins[1] - bins[0]
        counts, bins = np.histogram(data, bins=50, density=False)
        params_g = (np.max(counts), np.std(data), 0, 0)
        bounds = [(0, 2*params_g[0]), (0.0001, 5*params_g[1]), (None, None), (0, 2*params_g[0])]
        if fit_type == 'gauss':
            objective = self.error
        elif fit_type == 'lorentz':
            objective = self.error_lorentz
        elif fit_type == 'quad':
            objective = self.error_quad
            bounds.append((None, None))
            params_g = (np.max(counts), np.std(data), 0, 0, 0)
        x = (bins[1:] + bins[:-1])/2
        y_data = counts#/width
        tol = 0.0000001
        print(params_g)
        res = scipy.optimize.minimize(objective, params_g, args=[y_data, x], bounds=bounds, tol=tol)
        end_err = objective(res.x, [y_data, x])
        print(f'{objective(params_g, [y_data, x]):.3f} -> {end_err:.3f}')
        return res, end_err

    def fit_model(self, x, params, fit_type='gauss'):
        if fit_type == 'gauss':
            objective = self.model
        elif fit_type == 'lorentz':
            objective = self.model_lorentz
        elif fit_type == 'quad':
            objective = self.model_quad
        return objective(x, params)

    def fit_error(self, x, params, fit_type='gauss'):
        if fit_type == 'gauss':
            objective = self.error
        elif fit_type == 'lorentz':
            objective = self.error_lorentz
        elif fit_type == 'quad':
            objective = self.error_quad
        return objective(x, params)

    def dof(self, fit_type='gauss'):
        if fit_type == 'gauss':
            return 4
        elif fit_type == 'lorentz':
            return 4
        elif fit_type == 'quad':
            return 5

    @staticmethod
    def model(x, params):
        return params[0]*np.exp(-x**2/(2*params[1]**2)) + params[2]*x + params[3]

    @staticmethod
    def model_lorentz(x, params):
        return params[0]*(params[1]**2/(x**2 + params[1]**2)) + params[2]*x + params[3]

    @staticmethod
    def model_quad(x, params):
        return params[0]*np.exp(-x**2/(2*params[1]**2)) + params[2]*x + params[3] + params[4]*x**2
    
    @staticmethod
    def gauss(x, params):
        return params[0]*np.exp(-x**2/(2*params[1]**2))
    
    def error(self, params, args):
        return np.sum((args[0] - self.model(args[1], params))**2)

    def error_lorentz(self, params, args):
        return np.sum((args[0] - self.model_lorentz(args[1], params))**2)

    def error_quad(self, params, args):
        return np.sum((args[0] - self.model_quad(args[1], params))**2)

    def annulus_correction(self, z, z_ref=0.8):
        return


    def _inverse_comoving_radial_dist(self, r):
        '''
        Inverse of cosmo.comoving_radial_distance since tabulatedDESI does not have it
        '''
        return np.interp(r, self.cosmo._comoving_radial_distance, self.cosmo._z, left=None, right=None)

    def run(self, clusfrac_override=None, skip_background=False, physical=True, rbins=15, angbins=18, fit=False, fit_type='gauss'):
        if physical:
            # order here currently doesn't matter but it might be a good change to do a binning of sperp for each rbin
            self.bin_r(rbins)
            self.bin_sperp(angbins)
            if fit:
                mistab = self.impute_physical_fit(clusfrac_override=clusfrac_override, skip_background=skip_background, fit_type=fit_type)
            else:
                mistab = self.impute_physical(clusfrac_override=clusfrac_override, skip_background=skip_background)
        else:
            self.bin_angular()
            self.bin_z()
            mistab = self.impute(clusfrac_override=clusfrac_override, skip_background=skip_background)
        return mistab