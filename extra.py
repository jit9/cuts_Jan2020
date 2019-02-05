"""This script consists the scripts that are implemented in moby2
but are not used in generating the pickle parameters that we use
"""

class GetSlowMode(Routine):
    def __init__(self, **params):
        self._params = params
        self._tod = params.get('tod', None)
        self._driftFilter = params.get('driftFilter', None)
        self._output_key = paramsg.et('output_key', None)

    def execute(self, store):
        tod = store.get(self._tod)

        # get the range of frequencies to work on
        n_l = 1
        n_h = nextregular(int(round(self._driftFilter/df))) + 1

        # extract the fourior modes
        lf_data = fdata[:,n_l:n_h]

        # calculate the mean fourior modes for live and dark 
        fcmL = lf_data[self.preLiveSel].mean(axis = 0)
        fcmD = lf_data[self.preDarkSel].mean(axis = 0)

        # get the common modes for both live and dark
        dsCM, dsCM_dt = get_time_domain_modes(fcmL,n_l, tod.nsamps, df)
        dsDCM, _ = get_time_domain_modes(fcmD,n_l, tod.nsamps, df)

        results = {
            "dsCM": dsCM,
            "dsDCM": dsDCM
        }
        store.set(self._output_key, results)

    
class Retrend(Routine):
    def __init__(self, **params):
        self._params = params


    def execute(self, store):
        # get the trend for the live detectors
        trL = numpy.array(trend).T[self.preLiveSel].mean(axis=0)
        trLt = trL[:,np.newaxis]

        # retrend the live detectors
        moby2.tod.retrend_tod(trLt, data = self.dsCM)
        
        # get the trend for the dark detectors
        trD = numpy.array(trend).T[self.preDarkSel].mean(axis=0)
        trDt = trD[:,np.newaxis]

        # retrend the dark detectors
        moby2.tod.retrend_tod(trDt, data = self.dsDCM)


class AnalyzeAtm(Routine):
    def __init__(self, **params):
        """This routine does the atomosphere 1/f analysis"""
        self._params = params
        self._fitPowerLaw = params.get("fitPowerLaw", False)

    def execute(self, store):
        if self._fitPowerLaw:
            # look at both the live and dark detectors
            sel = preLiveSel + preDarkSel

            # combine the rms from both live and dark detectors
            rms = self.crit["rmsLive"]["values"] + self.crit["rmsDark"]["values"]

            # fit an atmosphere model with it 
            powLaw, level, knee = fit_atm(fdata, sel, dt, df, rms, self.scan_freq,
                                          **par.get("atmFit",{}))
            results = {
                "atmPowLaw": powLaw,
                "atmLevel": level,
                "atmKnee": knee
            }
        self.store(self._output_key, results)


            
def fit_atm(fdata, sel, dt, df, noise, scanf, 
            fminA=0.2, fmaxA=3., fmaxT = 10., 
            its = 1, width = 0.005):
    """
    Fit a power law to the atmosphere signal in a range of frequencies.
    """
    scale = 2*dt**2*df
    delta = 0.7
    kneeM = fmaxA+delta
    ind_ini = int(fminA/df)
    ps = np.power(np.abs(fdata[sel,:int(fmaxT/df)]),2)*scale
    # Get scan harmonics
    n_harm = int(np.ceil(fmaxT/scanf))
    i_harm = np.array(np.round(np.arange(1,n_harm+1)*scanf/df), dtype=int)
    i_harmw = i_harm.copy()
    di = int(width/df)
    for i in xrange(di):
        i_harmw = np.hstack([i_harmw,i_harm-(i+1)])
        i_harmw = np.hstack([i_harmw,i_harm+(i+1)])
    i_harmw.sort()
    # Iterate range fit
    for it in xrange(its):
        fmax = kneeM-delta-fminA
        imin = int(fminA/df); imax = int(fmax/df)
        psr = ps[:,imin:imax]
        log_ps = np.log(psr)
        freq = np.arange(imin,imax)*df
        log_freq = np.log(freq)
        w = np.diff(log_freq)
        w = np.hstack([w,w[-1]])
        iharm = i_harmw[(i_harmw>imin)*(i_harmw<imax)]-imin
        s = np.ones(len(freq),dtype=bool)
        s[iharm] = False
        m,n = np.polyfit(log_freq[s], log_ps[:,s].T, 1, w=w[s])
        pl = np.power(freq[np.newaxis,s].repeat(ps.shape[0],0),
                      m[:,np.newaxis].repeat(len(freq[s]),1))*\
                      np.exp(n[:,np.newaxis].repeat(len(freq[s]),1))
        c = np.sum(psr[:,s]*pl,axis=1)/np.sum(pl*pl,axis=1)
        level = np.exp(n)*c
        knee = np.power(noise[sel]/level,1./m)
        kneeM = np.median(knee)
    mA = np.zeros(fdata.shape[0])
    levelA = np.zeros(fdata.shape[0])
    kneeA = np.zeros(fdata.shape[0])
    mA[sel] = m
    levelA[sel] = level
    kneeA[sel] = knee
    return mA, levelA, kneeA            
