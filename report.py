from todloop import Routine
import cPickle, h5py, os
import numpy as np


class Summarize(Routine):
    def __init__(self, **params):
        """This routine aims to summarize the parameters derived
        from all the other routines"""
        Routine.__init__(self)
        self.inputs = params.get('inputs', None)
        self.outputs = params.get('outputs', None)

    def execute(self, store):
        # initialize an empty dictionary to store results
        results = {}
        
        # retrieve all the calculated results
        for key in self.inputs.keys():
            results.update(store.get(self.inputs[key]))

        self.logger.info("Successfully processed: %s" % results.keys())
        store.set(self.outputs.get('report'), results)

        
class PrepareDataLabel(Routine):
    def __init__(self, **params):
        """Prepare an HDF5 data set that contains all relevant metedata
        and detector timeseries as a basis for various ML studies
        """
        Routine.__init__(self)
        self.inputs = params.get("inputs", None)
        self._pickle_file = params.get('pickle_file', None)
        self._output_file = params.get('output_file', None)
        self._group_name = params.get('group', None)
        self._downsample = params.get('downsample', 1)
        
    def initialize(self):
        # load pickle file
        self.logger.info("Loading %s..." % self._pickle_file)
        with open(self._pickle_file, "r") as f:
            self._pickle_data = cPickle.load(f)

        # create output h5 file if it doesn't exist
        if os.path.isfile(self._output_file):
            # file exists
            self.logger.info("Found %s, updating instead" % self._output_file)
            self._hf = h5py.File(self._output_file, 'a')
        else:
            # file doesn't exist
            self.logger.info("Creating %s..." % self._output_file)
            self._hf = h5py.File(self._output_file, 'w')
            
        try:
            self.logger.info("Creating group %s..." % self._group_name)            
            self._group = self._hf.create_group(self._group_name)
        except ValueError:
            self.logger.info("%s exists, update instead" % self._group_name)            
            self._group = self._hf[self._group_name]
            
        
    def execute(self, store):
        # retrieve tod
        tod = store.get(self.inputs.get('tod'))

        # retrieve the calculated statustics
        report = store.get(self.inputs.get('report'))
        keys = report.keys()
        
        # get relevant metadata for this tod from pickle file
        tod_name = self.get_name()
        pickle_id = self._pickle_data['name'].index(tod_name)

        # get tes mask
        tes_mask = list(np.where(tod.info.array_data['det_type'] == 'tes'))[0]

        # store each det timeseries in hdf5
        for tes_det in tes_mask:
            if self._downsample:
                data = tod.data[tes_det, ::self._downsample]
            else:
                data = tod.data[tes_det]
        
            # generate a unique detector id
            det_uid = '%d.%d' % (self.get_id(), tes_det)

            # test if the dataset already exists, if so update it instead
            try:
                dataset = self._group.create_dataset(det_uid, data=data)
            except RuntimeError:
                dataset = self._group[det_uid]
                dataset[:] = data

            # save report to h5 file
            for k in keys:
                dataset.attrs[k] = report[k][tes_det]
                
            # save label
            dataset.attrs['label'] = int(self._pickle_data['sel'][tes_det, pickle_id])

        self.logger.info("Data saved in %s" % self._output_file)
 
    def finalize(self):
        self._hf.close()
        
