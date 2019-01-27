import os
import numpy as np
import scipy.stats.mstats as ms

import moby2
from moby2.scripting import products
from moby2.analysis import hwp
from todloop import Routine


class CutSources(Routine):
    def __init__(self, **params):
        """A routine that cuts the point sources"""
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._tag_source = params.get('tag_source', None)
        self._source_list = params.get('source_list', None)
        self._no_noise = params.get('no_noise', True)
        self._pointing_par = params.get('pointing_par', None)
        self._mask_params = params.get('mask_params', {})
        self._shift_params = params.get('mask_shift_generator', None)
        self._depot_path = params.get('depot', None)

    def initialize(self):
        self._depot = moby2.util.Depot(self._depot_path)
        user_config = moby2.util.get_user_config()
        moby2.pointing.set_bulletin_A(params=user_config.get('bulletin_A_settings'))

    def execute(self, store):
        # retrieve tod
        tod = store.get(self._input_key)
        
        # check if source cut results exist
        sourceResult = os.path.exists(
            self._depot.get_full_path(
                moby2.TODCuts, tag=self._tag_source, tod=tod))
        
        # if cuts exist, load it now
        if sourceResult:
            self.logger.info("Loading time stream cuts (%s)" % self._tag_source)

            # load source cut
            source_cuts = self._depot.read_object(
                moby2.TODCuts, tag=self._tag_source, tod=tod)
            # fill the cuts in the TOD
            moby2.tod.fill_cuts(tod, source_cuts, no_noise=self._no_noise)

        # if source cut cannot be retrieved by tag_source, load it
        # through _source_list
        elif self._source_list is not None:
            self.logger.info("Finding new source cuts")

            # retrieve source list from the file given
            with open(self._source_list, 'r') as f:
                source_list = f.readlines()
                source_list = [(s.strip('\n'), 'source') for s in source_list]

            # supply focal plane information to tod
            tod.fplane = products.get_focal_plane(self._pointing_par, tod.info)

            # find sources that fall in the given TOD
            matched_sources = moby2.ephem.get_sources_in_patch(
                tod=tod, source_list=source_list)

            # check if shift is needed
            if self._shift_params is not None:
                # calculate pointing offset
                offset = products.get_pointing_offset(
                    self._shift_params, tod=tod, source_offset=True)

                # check if offset is calculated successfully, if not give
                # a zero offset
                if offset is None:
                    offset = (0., 0.)

                # calculate a map size
                if max(offset) > 20. / 60:
                    self._mask_params['map_size'] = max(offset) + 10. / 60
                self._mask_params['offset'] = offset

            self.logger.info("matched sources: %s" % matched_sources)

            # create a placeholder cut object to store our source cuts
            pos_cuts_sources = moby2.TODCuts.for_tod(tod, assign=False)

            # process source cut for each source
            for source in matched_sources:
                # compute the source cut associated with the source
                source_cut = moby2.tod.get_source_cuts(
                    tod, source[1], source[2], **self._mask_params)
                # merge the source cut to the total cuts
                pos_cuts_sources.merge_tod_cuts(source_cut)

            # fill the source cuts to the tod
            moby2.tod.fill_cuts(tod, pos_cuts_sources, no_noise=self._no_noise)

            # write to depot, copied from moby2, not needed here
            # depot.write_object(pos_cuts_sources,
            #                    tag=params.get('tag_source'),
            #                    force=True, tod=tod, make_dirs=True)

        # pass the processed tod back to data store
        store.set(self._output_key, tod)


class CutPlanets(Routine):
    def __init__(self, **params):
        """A routine that perform the planet cuts"""
        Routine.__init__(self)        
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._no_noise = params.get('no_noise', True)
        self._tag_planet = params.get('tag_planet', None)
        self._pointing_par = params.get('pointing_par', None)
        self._mask_params = params.get('mask_params', {})
        self._shift_params = params.get('mask_shift_generator', None)
        self._depot_path = params.get('depot', None)

    def initialize(self):
        self._depot = moby2.util.Depot(self._depot_path)
        
    def execute(self, store):
        tod = store.get(self._input_key)

        # check if planetCuts exist
        planetResult = os.path.exists(
            self._depot.get_full_path(
                moby2.TODCuts, tag=self._tag_planet, tod=tod))
        

        # if planetCuts exist load it into variable pos_cuts_planets
        if planetResult:
            self.logger.info("Loading time stream cuts (%s)" % self._tag_planet)
            pos_cuts_planets = self._depot.read_object(
                moby2.TODCuts, tag=self._tag_planet, tod=tod)

        # if planetCuts do not exist generate it on the run
        else:
            self.logger.info("Finding new planet cuts")
            if not hasattr(tod, 'fplane'):
                tod.fplane = products.get_focal_plane(self._pointing_par,
                                                      tod.info)
            # load planet sources
            matched_sources = moby2.ephem.get_sources_in_patch(
                tod=tod, source_list=None)

            # check if shift is needed
            if self._shift_params is not None:
                # calculate pointing offset
                offset = products.get_pointing_offset(
                    self._shift_params, tod=tod, source_offset=True)

                # check if offset is calculated successfully, if not give
                # a zero offset
                if offset is None:
                    offset = (0., 0.)

                # calculate a map size
                if max(offset) > 20. / 60:
                    self._mask_params['map_size'] = max(offset) + 10. / 60
                self._mask_params['offset'] = offset

            self.logger.info("matched sources: %s" % matched_sources)

            # a place holder cut object to store all planet cut
            pos_cuts_planets = moby2.TODCuts.for_tod(tod, assign=False)

            # process each planet source
            for source in matched_sources:

                # calculate planet cut
                planet_cut = moby2.tod.get_source_cuts(
                    tod, source[1], source[2], **self._mask_params)
                # merge it into the total cut
                pos_cuts_planets.merge_tod_cuts(planet_cut)

            # write planet cut to depot, copied from moby2, not needed
            # here depot.write_object(pos_cuts_planets,
            # tag=params.get('tag_planet'), force=True, tod=tod,
            # make_dirs=True)

        # fill planet cuts into tod
        moby2.tod.fill_cuts(tod, pos_cuts_planets, no_noise=self._no_noise)

        # pass the processed tod back to data store
        store.set(self._output_key, tod)


class RemoveSyncPickup(Routine):
    def __init__(self, **params):
        """This routine fit / removes synchronous pickup"""
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._remove_sync = params.get('remove_sync', False)
        self._force_sync = params.get('force_sync', False)
        self._tag_sync = params.get('tag_sync', None)
        self._depot_path = params.get('depot', None)

    def initialize(self):
        self._depot = moby2.util.Depot(self._depot_path)

    def execute(self, store):
        # retrieve tod
        tod = store.get(self._input_key)

        # Check for existing results, to set what operations must be
        # done/redone.
        sync_result = os.path.exists(
            self._depot.get_full_path(
                moby2.tod.Sync, tag=self._tag_sync, tod=tod))

        # determine if sync is needed
        skip_sync = not self._remove_sync or (not self._force_sync
                                              and sync_result)

        # obtain scan frequency
        scan_freq = moby2.tod.get_scan_info(tod).scan_freq

        if (self._remove_sync) and (scan_freq != 0):
            self.logger.info("Removing Sync")
            # check if sync can be skipped
            if skip_sync:
                self.logger.info("Using old sync")
                ss = self._depot.read_object(
                    moby2.tod.Sync, tag=self._tag_sync, tod=tod)
            # if not generate it on the go
            else:
                self.logger.info("Computing new sync")
                ss = moby2.tod.Sync(tod)
                ss.findOutliers()
                ss = ss.extend()

                # write sync object to disk
                # depot.write_object(ss, tag=self._tag_sync, tod=tod, make_dirs=True,
                #                    force=True)

            ss.removeAll()
            del ss

        # pass the processed tod back to data store
        store.set(self._output_key, tod)


class CutPartial(Routine):
    def __init__(self, **params):
        """A routine that performs the partial cuts"""
        Routine.__init__(self)        
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._tag_partial = params.get('tag_partial', None)
        self._force_partial = params.get('force_partial', False)
        self._glitchp = params.get('glitchp', {})
        self._include_mce = params.get('include_mce', True)
        self._depot_path = params.get('depot', None)
        self._no_noise = params.get('no_noise', True)

    def initialize(self):
        self._depot = moby2.util.Depot(self._depot_path)

    def execute(self, store):
        # retrieve tod
        tod = store.get(self._input_key)

        # check if partial results already exist
        partial_result = os.path.exists(
            self._depot.get_full_path(moby2.TODCuts,
                                      tag=self._tag_partial, tod=tod))
        # check if we need to skip creating partial cuts
        skip_partial = not self._force_partial and partial_result

        # if we want to skip creating partial cuts, load from depot
        if skip_partial:
            # Read existing result
            self.logger.info("Loading time stream cuts (%s)" % self._tag_partial)            
            cuts_partial = self._depot.read_object(
                moby2.TODCuts, tag=self._tag_partial, tod=tod)
        # otherwise generate partial cuts now
        else:
            self.logger.info('Generating partial cuts')

            # Generate and save new glitch cuts (note calbol may not be implemented...)
            cuts_partial = moby2.tod.get_glitch_cuts(
                tod=tod, params=self._glitchp)

        # check if we want to include mce_cuts
        if self._include_mce:
            # find mce cuts
            mce_cuts = moby2.tod.get_mce_cuts(tod)

            # merge it with the partial cuts
            cuts_partial.merge_tod_cuts(mce_cuts)

            # write to depot, not needed here
            # depot.write_object(cuts_partial,
            # tag=params.get('tag_partial'), tod=tod, make_dirs =
            # True, force=True)

        # fill the partial cuts in our tod
        moby2.tod.fill_cuts(
            tod, cuts_partial, extrapolate=False, no_noise=self._no_noise)

        # save the partial cuts in tod object for further processing
        tod.cuts = cuts_partial

        # pass the tod back to the store
        store.set(self._output_key, tod)


class SubstractHWP(Routine):
    def __init__(self, input_key, output_key, **params):
        """This routine substracts the A(chi) signal from HWP"""
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._hwp_par = params.get('hwp_par')
        self._depot_path = params.get('depot', None)

    def initialize(self):
        self._depot = moby2.util.Depot(self._depot_path)

    def execute(self, store):
        # retrieve tod
        tod = store.get(self._input_key)

        self.logger.info("Substract HWP signal")

        # retrieve hwp_modes object from depot
        hwp_modes = self._depot.read_object(
            hwp.HWPModes,
            tag=self._hwp_par['a_chi']['tag'],
            tod=tod,
            structure=self._hwp_par['a_chi']['structure'])

        # get hwp angles
        hwp_angles = moby2.scripting.products.get_hwp_angles(
            self._hwp_par['angles'], tod)

        # substracting the hwp sinal
        r = hwp_modes.get_reconstructor(hwp_angles * np.pi / 180)
        hwp_signal = r.get_achi()
        tod.data[hwp_modes.det_uid, :] -= hwp_signal

        # pass the tod to the data store
        store.set(self._output_key, tod)


class TransformTOD(Routine):
    def __init__(self, **params):
        """This routine transforms a series of tod data transformation
        such as downsampling, remove_mean and detrend"""
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._remove_mean = params.get('remove_mean', True)
        self._remove_median = params.get('remove_mediam', False)
        self._detrend = params.get('detrend', True)
        self._remove_filter_gain = params.get('remove_filter_gain', False)
        self._n_downsample = params.get('n_downsample', None)

    def execute(self, store):
        # retrieve tod
        tod = store.get(self._input_key)

        # remove mean or remove median
        if self._remove_mean:
            moby2.tod.remove_mean(tod)
        elif self._remove_median:
            moby2.tod.remove_median(tod)

        # detrend
        if self._detrend:
            moby2.tod.detrend_tod(tod)

        # remove filter gain
        if self._remove_filter_gain:
            moby2.tod.remove_filter_gain(tod)

        # downsampling
        if self._n_downsample is not None:
            tod = tod.copy(resample=2**self._n_downsample, resample_offset=1)
            self.logger.info("Downsampling done")

        store.set(self._output_key, tod)


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

        chunk_params = {
            'T': scan["T"] * ds,
            'pivot': scan["pivot"] * ds,
            'N': scan["N"]
        }
        self.logger.info(chunk_params)
        store.set(self._output_key, chunk_params)

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
        
        if self._channel is None or self._T_max is None or self._dT_max is None:
            pass
        else:
            thermometers = []
            for ch in channels:
                thermometer = tod.get_hk( ch, fix_gaps=True)
                if len(np.diff(thermometer).nonzero()[0]) > 0:
                    thermometers.append(thermometer)
            if len(thermometers) > 0:
                thermometers = numpy.array(thermometers)
                
                # Get thermometer statistics
                th_mean = moby2.tod.remove_mean(data=thermometers)
                th_trend = moby2.tod.detrend_tod(data=thermometers)
                Temp = th_mean[0]
                dTemp = th_trend[1][0] - th_trend[0][0]
                if (Temp > T_max) or (abs(dTemp) > dT_max):
                    temperatureCut = True        

        thermal_results = {
            'Temp': Temp,
            'dTemp': dTemp,
            'temperatureCut': temperatureCut
        }
        
        self.logger.info(thermal_results)
        store.set(self._output_key, thermal_results)


class GetDetectors(Routine):
    def __init__(self, **params):
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._fullRMSlim = params.get('fullRMSlim', 1e8)
        self._source = params.get('source', None)
        self._filename = params.get('filename', None)
        self._live = params.get('live', None)
        self._dark = params.get('dark', None)
        self._exclude = params.get('exclude', None)
        self._noExclude = params.get('noExclude', False)

    def execute(self, store):
        # get tod
        tod = store.get(self._input_key)

        # get all detector lists
        # directly copied from loic, not sure why copy is needed here
        # i think it might be just a safety precaution
        dets = tod.info.det_uid().copy()

        # retrieve predefined exclude, dark live detector lists from file
        _exclude, _dark, _live = self.get_detector_params()

        # create mask based on the provided lists

        # look at the exclude detector list
        # if the given list is based on individual source
        if _exclude.has_key('det_uid'):
            exclude = tod.info.array_data.select_inner({
                'det_uid': _exclude['det_uid']
            }, mask = True, det_uid = dets)
            
        # if the given list is based on matrix source
        else:
            exclude = tod.info.array_data.select_inner({
                'row': exclude['rows'],
                'col': exclude['cols']
            }, mask = True, det_uid = dets)

        # noExclude parameter specifies if we want to remove
        # the detectors that are cutted beforehand, if it's
        # set as false, then the existing cuts are retrieved
        # through the get_cut method and marked as exclude
        # on top of what's given
        if not(self._noExclude):
            exclude[list(tod.cuts.get_cut())] = True

        # look at dark detector list
        # if the list is based on individual source
        if _dark.has_key('det_uid'):
            dark_candidates = tod.info.array_data.select_inner({
                'det_uid': _dark['det_uid']
            }, mask = True, det_uid = dets)
            
        # if the list is based on matrix source
        else:
            dark_candidates = tod.info.array_data.select_inner({
                'row': _dark['rows'],
                'col': _dark['cols']
            }, mask = True, det_uid = dets)

        # look at the live detector lists
        # if the list is based on individual source
        if _live.has_key('det_uid'):
            live_candidates = tod.info.array_data.select_inner({
                'det_uid': _live['det_uid']
            }, mask = True, det_uid = dets)
            
        # if the list is based on matrix source
        else:
            live_candidates = tod.info.array_data.select_inner({
                'row': _live['rows'],
                'col': _live['cols']
            }, mask = True, det_uid = self.dets)
    
        # filter zero detectors
        # mark zero detectors as 1, otherwise 0 
        self.logger.info('Finding zero detectors')
        zero_sel = ~tod.data[:,::100].any(axis=1)

        # filter detectors with too large rms
        # mark good detectors as 1 and bad (large rms) as 0
        full_rms_sel = np.std(tod.data, axis=1) < self._fullRMSlim

        # exclude zero detectors and noisy detectors
        live = live_candidates * ~zero_sel * full_rms_sel
        dark = orig_candidates * ~zero_sel * full_rms_sel

        # export the detectors lists
        detectors = {
            'live_candidates': live_candidates,
            'dark_candidates': dark_candidates,
            'live_final': live,
            'dark_final': dark
        }
        self.logger.info(detectors)
        store.set(self._output_key, detectors)

    def get_detector_params(self):
        source = self._source

        # if given detector lists are specified in terms of row / col matrix
        # so far I haven't seen it being used
        if source == "matrix":
            mat = np.loadtxt(self._filename, dtype = str)
            r,c = np.where(mat == "d")  # Dark detectors
            dark = {}.fromkeys(("rows", "cols"))
            dark["rows"] = r; dark["cols"] = c
            r,c = np.where(mat == "b")  # Bad detectors
            exclude = {}.fromkeys(("rows", "cols"))
            exclude["rows"] = r; exclude["cols"] = c
            r, c = np.where(mat == "s")  # Stable detectors
            r2, c2 = np.where(mat == "c")  # Candidate detectors
            liveCandidates = {}.fromkeys(("rows", "cols"))
            liveCandidates["rows"] = np.hstack([r,r2])
            liveCandidates["cols"] = np.hstack([c,c2])

        # if given detector lists are specified in term of lists
        elif source == "individual":
            # load excludable detector list
            exclude = moby2.util.MobyDict.from_file(self._exclude))
            
            # load dark detector list
            dark = moby2.util.MobyDict.from_file(self._dark)
            
            # load live detectors
            liveCandidates = moby2.util.MobyDict.from_file(self._live)
            
        else:
            raise "Unknown detector params source"
        
        return exclude, dark, liveCandidates        


class CalibrateTOD(Routine):
    def __init__(self, **params):
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._dets_key = params.get('dets_key', None)
        self._output_key = params.get('output_key', None)
        self._output_calData = params.get('output_calData', "calData")
        self._flatfield = params.get('flatfield', None)
        self._forceNoResp = params.get('forceNoResp', None)
        self._config = params.get('config', None)
        self._calibrateTOD = params.get('calibrateTOD', True)

    def execute(self, store):
        tod = store.get(self._input_key)

        #####################################################
        # get responsivities and flatfield for calibration  #
        #####################################################
        
        # get responsivity
        resp = products.get_calibration(self._config, tod.info)

        # select only responsive detectors
        respSel = (resp.cal != 0.0)
        
        # get flatfield and a selection mask
        flatfield_object = moby2.detectors.RelCal.from_dict(self._flatfield)
        ffSel, ff = flatfield_object.get_property('cal',
                                                  det_uid=tod.info.dets,
                                                  default=1.)
        
        # get stable detectors
        _, stable = flatfield_object.get_property('stable',
                                                  det_uid=tod.info.dets,
                                                  default=False)

        # check if we want to fill default responsivity
        if self._forceNoResp:
            # fill the default with median of stable detectors
            rm = np.median(resp.cal[stable*respSel])
            resp.cal[~respSel] = rm

        # get the RMS for flatfield calibration if it exists
        # otherwise fill with 0
        if flatfield_object.calRMS is not None:
            _, ffRMS = flatfield_object.get_property('calRMS',
                                                     det_uid=tod.info.dets,
                                                     default=1.)
        else:
            ffRMS = np.zeros_like(tod.info.dets)

        # summarize all the calibration data into a dictionary 
        calData = {
            "resp": resp.cal,
            "respSel": respSel,            
            "ff": ff,
            "ffRMS": ffRMS,
            "ffSel": ffSel,
            "stable": stable,
            "cal": resp.cal*ff,
            "calSel": ffSel*respSel,
            "calibrated": False
        }

        
        ########################################################
        # calibrate TOD to pW using responsivity and flatfield #
        ########################################################
        
        cal = calData['cal']
        
        # apply to all except for original dark detectors        
        orig_dark = store.get(self._dets_key)['dark_candidates']
        s = ~orig_dark
        
        if not self._calibratedTOD:
            moby2.libactpol.apply_calibration(tod.data,
                                              s.nonzero()[0].astype('int32'),
                                              cal[s].astype('float32'))
            # mark tod as calibrated
            calData['calibrated'] = True
            
        # report error if calibration is unsuccessful
        if not(np.any(self.calData["calSel"])): 
            self.logger.error('moby', 0, "ERROR: no calibration for this TOD") 
            return 1

        # save to data store
        store.set(self._output_calData, calData)
        store.set(self._output_key, tod)

