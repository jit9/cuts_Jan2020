import os
import numpy as np

from todloop import Routine

import moby2
from moby2.scripting import products
from moby2.analysis import hwp

class CutSources(Routine):
    def __init__(self, input_key, output_key, **params):
        """A routine that cuts the point sources"""
        Routine.__init__(self)
        self._input_key = input_key
        self._output_key = output_key
        self._tag_source = params.get('tag_source', None)
        self._source_list = params.get('source_list', None)
        self._no_noise = params.get('no_noise', True)
        self._pointing_par = params.get('pointing_par', None)
        self._mask_params = params.get('mask_params', {})
        self._shift_params = params.get('mask_shift_generator', None)
        self._depot_path = params.get('depot', None)

    def initialize(self):
        self._depot = moby2.util.Depot(self._depot_path)

    def execute(self, store):
        # retrieve tod
        tod = store.get(self._input_key)

        # check if source cut results exist
        sourceResult = os.path.exists(
            self._depot.get_full_path(moby2.TODCuts,
                                      tag=self._tag_source, tod=tod) ) 
        # if cuts exist, load it now
        if sourceResult:
            self.logger.trace(0, "Loading time stream cuts (%s)" %
                              self._tag_source)

            # load source cut
            source_cuts = self._depot.read_object(moby2.TODCuts,
                                                  tag=self._tag_source,
                                                  tod=tod)
            # fill the cuts in the TOD
            moby2.tod.fill_cuts(tod, source_cuts,
                                no_noise=self._no_noise)

        # if source cut cannot be retrieved by tag_source, load it
        # through _source_list
        elif self._source_list is not None:
            self.logger.trace(0, "Finding new source cuts")

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
                offset = products.get_pointing_offset(self._shift_params,
                                                      tod=tod,
                                                      source_offset=True)

                # check if offset is calculated successfully, if not give
                # a zero offset
                if offset is None:
                    offset = (0.,0.)

                # calculate a map size
                if max(offset) > 20./60:
                    self._mask_params['map_size'] = max(offset) + 10./60
                self._mask_params['offset'] = offset

            self.logger.info("matched sources: %s" % matched_sources)

            # create a placeholder cut object to store our source cuts
            pos_cuts_sources = moby2.TODCuts.for_tod(tod, assign=False)

            # process source cut for each source
            for source in matched_sources:
                # compute the source cut associated with the source
                source_cut = moby2.tod.get_source_cuts( tod, source[1],
                                                        source[2],
                                                        **self._mask_params)
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
    def __init__(self, input_key, output_key, **params):
        """A routine that perform the planet cuts"""
        self._input_key = input_key
        self._output_key = output_key
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
        planetResult = os.path.exists( self._depot.get_full_path(moby2.TODCuts,
                                                                 tag=self._tag_planet,
                                                                 tod=tod) ) 
        # if planetCuts exist load it into variable pos_cuts_planets
        if planetResult:
            self.logger.trace(0, "Loading time stream cuts (%s)" % self._tag_planet)
            pos_cuts_planets = self._depot.read_object(moby2.TODCuts,
                                                       tag=self._tag_planet, tod=tod)

        # if planetCuts do not exist generate it on the run
        else:
            logger.trace(0, "Finding new planet cuts")
            if not hasattr(tod,'fplane'):
                tod.fplane = products.get_focal_plane(self._pointing_par, tod.info)

            # load planet sources
            matched_sources = moby2.ephem.get_sources_in_patch( tod=tod,
                                                                source_list=None)

            # check if shift is needed
            if self._shift_params is not None:
                # calculate pointing offset
                offset = products.get_pointing_offset(self._shift_params,
                                                      tod=tod,
                                                      source_offset=True)

                # check if offset is calculated successfully, if not give
                # a zero offset
                if offset is None:
                    offset = (0.,0.)

                # calculate a map size
                if max(offset) > 20./60:
                    self._mask_params['map_size'] = max(offset) + 10./60
                self._mask_params['offset'] = offset

            self.logger.info("matched sources: %s" % matched_sources)

            # a place holder cut object to store all planet cut
            pos_cuts_planets = moby2.TODCuts.for_tod(tod, assign=False)

            # process each planet source
            for source in matched_sources:

                # calculate planet cut
                planet_cut = moby2.tod.get_source_cuts( tod, source[1],
                                                        source[2],
                                                        **self._mask_params)
                # merge it into the total cut
                pos_cuts_planets.merge_tod_cuts(planet_cut)

            # write planet cut to depot, copied from moby2, not needed here
            # depot.write_object(pos_cuts_planets, tag=params.get('tag_planet'),
            #                    force=True, tod=tod, make_dirs=True)

        # fill planet cuts into tod
        moby2.tod.fill_cuts(tod, pos_cuts_planets, no_noise=self._no_noise)

        # pass the processed tod back to data store
        store.set(self._output_key, tod)


class RemoveSyncPickup(Routine):
    def __init__(self, input_key, output_key, **params):
        """This routine fit / removes synchronous pickup"""
        Routine.__init__(self)
        self._input_key = input_key
        self._output_key = output_key
        self._remove_sync = params.get('remove_sync', False)
        self._force_sync = params.get('force_sync', False)
        self._tag_sync = params.get('tag_sync', None)
        self._depot_path = params.get('depot', None)
        
    def initialize(self):
        self._depot = moby2.util.Depot(self._depot_path)

    def execute(self, store):
        # retrieve tod
        tod = store.get(self._input_key)

        # Check for existing results, to set what operations must be done/redone.
        sync_result = os.path.exists( self._depot.get_full_path(moby2.tod.Sync,
                                                                tag=self._tag_sync,
                                                                tod=tod))

        # determine if sync is needed
        skip_sync = not self._remove_sync or ( not self._force_sync and
                                               sync_result )

        # obtain scan frequency
        scan_freq = moby2.tod.get_scan_info(tod).scan_freq

        self.logger.trace(0, "Removing Sync")

        # check if sync can be skipped
        if skip_sync:
            self.logger.trace(2, "Using old sync")
            ss = self._depot.read_object(moby2.tod.Sync,
                                         tag=self._tag_sync, tod=tod)
        # if not generate it on the go
        else:
            self.logger.trace(2, "Computing new sync")
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
    def __init__(self, input_key, output_key, **params):
        """A routine that performs the partial cuts"""
        self._input_key = input_key
        self._output_key = output_key
        self._tag_partial = params.get('tag_partial', None)
        self._force_partial = params.get('force_partial', False)
        self._skip_sync = params.get('skip_sync', False)
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
        partial_result = os.path.exists( depot.get_full_path(moby2.TODCuts,
                                                             tag=self._tag_partial,
                                                             tod=tod))
        # check if we need to skip creating partial cuts
        skip_partial = self._skip_sync and not self._force_partial and partial_result

        # if we want to skip creating partial cuts, load from depot
        if skip_partial:
            # Read existing result
            cuts_partial = self._depot.read_object(moby2.TODCuts,
                                                   tag=self._tag_partial,
                                                   tod=tod)
        # otherwise generate partial cuts now
        else:
            self.logger.info('Generating partial cuts')
            
            # Generate and save new glitch cuts (note calbol may not be implemented...)
            cuts_partial = moby2.tod.get_glitch_cuts(tod=tod, params=self._glitchp)
            
        
        # check if we want to include mce_cuts
        if self._include_mce:
            # find mce cuts
            mce_cuts = moby2.tod.get_mce_cuts(tod)

            # merge it with the partial cuts
            cuts_partial.merge_tod_cuts(mce_cuts)

            # write to depot, not needed here
            # depot.write_object(cuts_partial, tag=params.get('tag_partial'), tod=tod,
            #                    make_dirs = True,
            #                    force=True)

        # fill the partial cuts in our tod
        moby2.tod.fill_cuts(tod, cuts_partial, extrapolate = False, no_noise=self._no_noise)

        # save the partial cuts in tod object for further processing
        tod.cuts = cuts_partial

        # pass the tod back to the store
        store.set(self._output_key, tod)


class SubstractHWP(Routine):
    def __init__(self, input_key, output_key, **params):
        Routine.__init__(self)
        self._input_key = input_key
        self._output_key = output_key
        self._hwp_par = params.get('hwp_par')
        self._depot_path = params.get('depot', None)

    def initialize(self):
        self._depot = moby2.util.Depot(self._depot_path)

    def execute(self, store):
        # retrieve tod
        tod = store.get(self._input_key)

        self.logger.trace(0, "Substract HWP signal")

        # retrieve hwp_modes object from depot
        hwp_modes = self._depot.read_object(
            hwp.HWPModes,
            tag = self._hwp_par['a_chi']['tag'],
            tod=tod,
            structure= self._hwp_par['a_chi']['structure'] )

        # get hwp angles
        hwp_angles = moby2.scripting.products.get_hwp_angles(self._hwp_par['angles'], tod)

        # substracting the hwp sinal
        r = hwp_modes.get_reconstructor(hwp_angles*np.pi/180)
        hwp_signal = r.get_achi()
        tod.data[hwp_modes.det_uid,:]-= hwp_signal

        # pass the tod to the data store
        store.set(self._output_key, tod)
