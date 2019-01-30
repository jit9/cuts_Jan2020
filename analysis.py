import numpy as np
from numpy import ma
import scipy.stats.mstats as ms
from scipy.cluster.vq import kmeans2

import moby2
from todloop import Routine

from utils import *


class AnalyzeScan(Routine):
    def __init__(self, **params):
        """This routine analyzes the scan pattern"""
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._scan_params = params.get('scan_param', {})

    def execute(self, store):
        # load tod
        tod = store.get(self._input_key)

        sample_time = (tod.ctime[-1] - tod.ctime[0]) / (tod.ctime.shape[0]-1)
        # analyze scan and save result into a dictionary
        scan = self.analyze_scan(
            np.unwrap(tod.az), sample_time,
            **self._scan_params)

        # get scan frequency
        scan_freq = scan["scan_freq"]

        # get downsample level
        ds = tod.info.downsample_level

        # summary of scan parameters
        scan_params = {
            'T': scan["T"] * ds,
            'pivot': scan["pivot"] * ds,
            'N': scan["N"],
            'scan_freq': scan_freq
        }
        
        self.logger.info(scan_params)
        store.set(self._output_key, scan_params)

    def analyze_scan(self, az, dt=0.002508, N=50, vlim=0.01, qlim=0.01):
        """Find scan parameters and cuts"""

        # Find no motion

        # compute the 1% and 99% quantiles
        lo, hi = ms.mquantiles(az, (qlim, 1 - qlim))

        # compute the scan speed

        # compute the az steps
        daz = np.diff(az)

        # form the scan speed array
        # note that:
        # daz[0] = az_1 - az_0
        # daz[1] = az_2 - az_1
        # 2*daz[0] - daz[1] = 2 az_1 - 2 az_0 - az_2 + az_1
        #                   = 3 az_1 - 2 az_0 - az_2
        # this is just an estimate of the scan speed at t=0
        v_scan = np.r_[2 * daz[0] - daz[1], daz]

        # smooth the scan speed vector with a simple moving average of
        # length N, the last indexing is to ensure that the size of
        # the array is the same as v_scan
        v_smooth = np.convolve(v_scan,
                               np.ones(N) / N)[(N - 1) / 2:-(N - 1) / 2]

        # estimate the speed using the median of the scan speeds
        speed = np.median(abs(v_smooth))

        # identify when the scanwhen the speed is either too fast
        # or two slow. The minimum speed requirement is specified
        # as a fraction of the median speed using ~vlim~
        stop = abs(v_scan) < vlim * speed
        pick = abs(v_smooth) > 2 * speed

        # exit now in case of a stare TOD
        # i doubt that this is ever going to occur
        if all(stop):
            scan = {
                "az_max": hi,
                "az_min": lo,
                "az_speed": speed,
                "scan_freq": 0.0,
                "az_cuts": None,
                "T": len(az),
                "pivot": 0,
                "N": 1
            }
            return scan

        # Find sections with no scanning by identifying
        # when the speed is below threshold and when the
        # scan range is an outlier
        noscan = stop * (az > lo) * (az < hi)

        # Get scan frequency
        # first calculate the fourior transform
        faz = np.fft.rfft(az - az.mean())
        # identify the highest frequency which corresponds to the
        # scan frequency
        fscan = np.where(abs(faz) == abs(faz).max())[0][0] / dt / len(az)

        # Find turnarounds
        az_min = np.median(az[stop * (az < lo)])
        az_max = np.median(az[stop * (az > hi)])
        # i don't understand this part
        td = abs(lo - az_min) + abs(az_max - hi)

        # Find scan period parameters
        T_scan = int(1. / fscan / dt)  # number of samples in scan period

        # this is kind of arbitrary
        T_ex = int(1.2 * T_scan)

        # find non-stopping part of scan
        onescan = az[~noscan][:T_ex]
        az_min0 = np.min(onescan[stop[~noscan][:T_ex] * (onescan < lo)
                                 * (onescan > lo - td)])
        az_max0 = np.max(onescan[stop[~noscan][:T_ex] * (onescan > hi)
                                 * (onescan < hi + td)])
        imin0 = np.where(az == az_min0)[0][0]
        imax0 = np.where(az == az_max0)[0][0]
        pivot = np.min([imin0, imax0])  # Index of first scan minima or maxima
        N_scan = (len(az) - pivot) / T_scan  # Number of complete scan periods

        # Find cuts
        if hi - lo < 1. * np.pi / 180:
            flag = np.ones_like(az, dtype=bool)
        else:
            flag = (pick + stop) * (az > lo) * (az < hi) + (az < lo - td) + (
                az > hi + td)

        c_vect = moby2.tod.cuts.CutsVector.from_mask(flag).get_buffered(100)

        # return scan parameters
        scan = {
            "az_max": az_max,
            "az_min": az_min,
            "az_speed": speed,
            "scan_freq": fscan,
            "az_cuts": c_vect,
            "T": T_scan,
            "pivot": pivot,
            "N": N_scan
        }
        return scan


class AnalyzeTemperature(Routine):
    def __init__(self, **params):
        """This routine will analyze the temperature of the TOD"""
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._channel = params.get('channel', None)
        self._T_max = params.get('T_max', False)
        self._dT_max = params.get('dT_max', None)

    def execute(self, store):
        """
        @brief Measure the mean temperature and thermal drift, and
               suggest a thermalCut
        @return mean temperature, thermal drift and thermal cut flag
        """
        tod = store.get(self._input_key)

        Temp = None
        dTemp = None
        temperatureCut = False
        
        if self._channel is None or self._T_max is None \
           or self._dT_max is None:
            pass
        else:
            thermometers = []
            for ch in self._channel:
                thermometer = tod.get_hk(ch, fix_gaps=True)
                if len(np.diff(thermometer).nonzero()[0]) > 0:
                    thermometers.append(thermometer)
            if len(thermometers) > 0:
                thermometers = np.array(thermometers)
                
                # Get thermometer statistics
                th_mean = moby2.tod.remove_mean(data=thermometers)
                th_trend = moby2.tod.detrend_tod(data=thermometers)
                Temp = th_mean[0]
                dTemp = th_trend[1][0] - th_trend[0][0]
                if (Temp > self._T_max) or (abs(dTemp) > self._dT_max):
                    temperatureCut = True
                    
        thermal_results = {
            'Temp': Temp,
            'dTemp': dTemp,
            'temperatureCut': temperatureCut
        }
        
        self.logger.info(thermal_results)
        store.set(self._output_key, thermal_results)


class FouriorTransform(Routine):
    def __init__(self, **params):
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._fft_data = params.get('fft_data', None)

    def execute(self, store):
        tod = store.get(self._input_key)

        # first de-trend tod 
        trend = moby2.tod.detrend_tod(tod)

        # find the next regular, this is to make fft faster
        nf = nextregular(tod.nsamps)
        fdata = np.fft.rfft(tod.data, nf)

        # time and freq units
        dt = (tod.ctime[-1]-tod.ctime[0])/(tod.nsamps-1)
        df = 1./(dt*nf)

        # summarize fft data
        fft_data = {
            'trend': trend,
            'fdata': fdata,
            'dt': dt,
            'df': df,
            'nf': nf
        }

        # store data into data store
        store.set(self._output_key, tod)
        store.set(self._fft_data, fft_data)


class AnalyzeDarkLF(Routine):
    def __init__(self, **params):
        self._dets = params.get('dets', None)
        self._fft_data = params.get('fft_data', None)
        self._tod = params.get('tod', None)
        self._output_key = params.get('output_key', None)
        self._scan = params.get('scan', None)
        self._freqRange = params.get('freqRange', None)
        self._params = params

    def execute(self, store):
        # retrieved relevant data from data store
        tod = store.get(self._tod)
        fft_data = store.get(self._fft_data)
        fdata = fft_data['fdata']
        df = fft_data['df']
        sel = store.get(self._dets)['dark_final']
        scan_freq = store.get(self._scan)['scan_freq']

        # get the frequency band parameters
        frange = self._freqRange
        fmin = frange.get("fmin", 0.017)
        fshift = frange.get("fshift", 0.009)
        band = frange.get("band", 0.070)
        Nwin = frange.get("Nwin", 1)

        psel = []
        corr = []
        gain = []
        norm = []

        minFreqElem = 16
        # loop over freq windows
        for i in xrange(Nwin):
            # find upper / lower bounds' corresponding index in freq
            # lower bound: fmin + [ fshifts ]
            # upper bound: fmin + band  + [ fshifts ]
            n_l = int(round((fmin + i*fshift)/df))
            n_h = int(round((fmin + i*fshift + band)/df))

            # if there are too few elements then add a few more to
            # have exactly the minimum required
            if (n_h - n_l) < minFreqElem:
                n_h = n_l + minFreqElem

            # perform low frequency analysis
            r = self.lowFreqAnal(fdata, sel, [n_l, n_h], df,
                                 tod.nsamps, scan_freq)

            # append the results to the relevant lists
            psel.append(r["preSel"])
            corr.append(r["corr"])
            gain.append(np.abs(r["gain"]))
            norm.append(r["norm"])

        # count pre-selected
        psel = np.sum(psel, axis=0)
        Nmax = psel.max()

        # number of pre-selected above median mask this marks a good
        # selection of detectors
        psel50 = psel >= Nmax/2.

        # Normalize gain by the average gain of a good selection of
        # detectors, here this selection is given by psel50 AND presel
        for g, s in zip(gain, psel):
            g /= np.mean(g[psel50*s])
        gain = np.array(gain)
        # give a default gain of 0 for invalid data
        gain[np.isnan(gain)] = 0.

        # use mean as representative values for gain
        mgain = ma.MaskedArray(gain, ~np.array(psel))
        mgain_mean = mgain.mean(axis=0)

        # use max as representative values for corr
        mcorr = ma.MaskedArray(corr, ~np.array(psel))
        mcorr_max = mcorr.max(axis=0)

        # use mean as representative values for norm
        mnorm = ma.MaskedArray(norm, ~np.array(psel))
        mnorm_mean = mnorm.mean(axis=0)
        
        # export the values
        results = {}
        
        results["corrDark"] = mcorr_max.data,
        results["gainDark"] = mgain_mean.data
        results["normDark"] = mnorm_mean.data
        results["darkSel"] = psel50.copy()  # not sure why copy is needed
        store.set(self._output_key, results)
        
    def lowFreqAnal(self, fdata, sel, frange, df, nsamps, scan_freq):
        """Find correlations and gains to the main common mode over a
        frequency range
        """
        # get relevant low freq data in the detectors selected
        lf_data = fdata[sel, frange[0]:frange[1]]
        ndet = len(sel)
        res = {}

        # Apply sine^2 taper to data
        if self._params.get("useTaper", False):
            taper = get_sine2_taper(frange, edge_factor = 6)
            lf_data *= np.repeat([taper],len(lf_data),axis=0)

        # Scan frequency rejection
        if self._params.get("cancelSync",False) and (scan_freq/df > 7):
            i_harm = get_iharm(frange, df, scan_freq,
                               wide=self._params.get("wide",True))
            lf_data[:, i_harm] = 0.0

        # Get correlation matrix
        c = np.dot(lf_data, lf_data.T.conjugate())
        a = np.linalg.norm(lf_data, axis=1)
        aa = np.outer(a,a)
        aa[aa==0.] = 1.
        cc = c/aa

        # Get Norm
        ppar = self._params.get("presel",{})
        norm = np.zeros(ndet,dtype=float)
        fnorm = np.sqrt(np.abs(np.diag(c)))
        norm[sel] = fnorm*np.sqrt(2./nsamps)
        nnorm = norm/np.sqrt(nsamps)

        # get a range of valid norm values 
        nlim = ppar.get("normLimit",[0.,1e15])
        if np.ndim(nlim) == 0:
            nlim = [0, nlim]
        normSel = (nnorm > nlim[0])*(nnorm < nlim[1])

        # If sigmaSep is specified, check if norms are divided in 2 groups,
        # and use the higher norms
        sigs = ppar.get("sigmaSep", None)
        if sigs is not None:
            # use k-means clustering to split the norm into two groups
            # cent refers to the center of each cluster and
            # lab refers to the label of each data (1: cluster 1; 2: cluster 2)
            cent, lab = kmeans2(nnorm[normSel], 2)

            # if all groups are larger than 20% of all data
            frac = 0.2
            if lab.sum() > len(lab)*frac and lab.sum() < len(lab)*(1-frac):
                # sort the norm value for both of the group
                c0 = np.sort(nnorm[normSel][lab==0])
                c1 = np.sort(nnorm[normSel][lab==1])

                # find the medium
                # not sure why this is needed versus just calling the medium
                mc0 = c0[len(c0)/2]
                mc1 = c1[len(c1)/2]

                # estimating the std using the fact that the std is
                # 0.741 times the interquantile range (1st and 3rd)
                # not sure why this is needed versus just calling the std
                sc0 = 0.741*(c0[(3*len(c0))/4] - c0[len(c0)/4])
                sc1 = 0.741*(c1[(3*len(c1))/4] - c1[len(c1)/4])

                # This calculation doesn't make sense to me, and it's probably
                # incorrect consider when mc0>mc1 this formula means nothing
                # [original formula]
                # sep =  (mc0 + sigs*sc0 - (mc1 - sigs*sc1))*np.sign(mc1-mc0)
                # [new formula]
                sigs_data = np.abs(mc0-mc1)/(sc0+sc1)
                if sigs_data > sigs:
                    # use the higher norm group
                    if mc1 > mc0:
                        normSel[normSel] *= (lab==1)
                    else:
                        normSel[normSel] *= (lab==0)
                        
            # in all other cases except when there is only one group
            # use the larger group, the other group is treated as an
            # outlier
            elif lab.sum() > 0:
                if lab.sum() > len(lab)/2:
                    normSel[normSel] *= (lab==1)
                else:
                    normSel[normSel] *= (lab==0)

        # check which preselection is specified
        presel_method = ppar.get("method", "median")
        if presel_method is "median":
            sl = presel_by_median(cc, sel=normSel[sel], **presel_params)
            res["groups"] = None
            
        elif presel_method is "groups":
            G, ind, ld, smap = group_detectors(cc, sel=normSel[sel], **presel_params)
            sl = np.zeros(cc.shape[1], dtype=bool)
            sl[ld] = True
            res["groups"] = {
                "G": G,
                "ind": ind,
                "ld": ld,
                "smap": smap
            }
        else:
            raise "ERROR: Unknown preselection method"

        # The number of sels are just overwhelmingly confusing
        # To clarify for myself,
        # - normSel: selects the detectors with good norm
        # - sel: the initial selection of detectors specified
        #        for this case it is selection of dark detectors
        # - sl: is the preselected detectors from the median
        #       or group methods
        # Here it's trying to apply the preselection to the
        # dark selection
        preSel = sel.copy()
        preSel[sel] = sl
        
        # Get Correlations
        u, s, v = np.linalg.svd(lf_data[sl], full_matrices=False )
        corr = np.zeros(ndet)
        if par.get("doubleMode", False):
            corr[preSel] = np.sqrt(abs(u[:,0]*s[0])**2+abs(u[:,1]*s[1])**2)/fnorm[sl]
        else:
            corr[preSel] = np.abs(u[:,0])*s[0]/fnorm[sl]

        # Get Gains
        # data = CM * gain
        gain = np.zeros(ndet, dtype=complex)
        gain[preSel] = np.abs(u[:, 0])
        res.update({"preSel": preSel, "corr": corr, "gain": gain, "norm": norm, 
                    "cc": cc, "normSel": normSel})
        
        return res


        
