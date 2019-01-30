import numpy as np

import moby2
from moby2.scripting import products
from todloop import Routine

from utils import *


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
        dets = tod.info.det_uid.copy()

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
        dark = dark_candidates * ~zero_sel * full_rms_sel

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
            exclude = moby2.util.MobyDict.from_file(self._exclude)
            
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
                                                  det_uid=tod.info.det_uid,
                                                  default=1.)
        
        # get stable detectors
        _, stable = flatfield_object.get_property('stable',
                                                  det_uid=tod.info.det_uid,
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
                                                     det_uid=tod.info.det_uid,
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
            "calibrated": False,
            "flatfield_object": flatfield_object
        }

        
        ########################################################
        # calibrate TOD to pW using responsivity and flatfield #
        ########################################################
        
        cal = calData['cal']
        
        # apply to all except for original dark detectors        
        orig_dark = store.get(self._dets_key)['dark_candidates']
        s = ~orig_dark
        
        if not self._calibrateTOD:
            moby2.libactpol.apply_calibration(tod.data,
                                              s.nonzero()[0].astype('int32'),
                                              cal[s].astype('float32'))
            # mark tod as calibrated
            calData['calibrated'] = True
            
        # report error if calibration is unsuccessful
        if not(np.any(calData["calSel"])): 
            self.logger.error('moby', 0, "ERROR: no calibration for this TOD") 
            return 1

        # save to data store
        store.set(self._output_calData, calData)
        store.set(self._output_key, tod)
        
