from todloop import Routine

class Summarize(Routine):
    def __init__(self, **params):
        """This routine aims to summarize the parameters derived
        from all the other routines and save it into a file for
        furthur analysis"""
        self.inputs = params.get('inputs', None)
        self.outpath = params.get('outpath', None)

    def execute(self, store):
        # retrieve all the calculated results
        lf_dark = store.get(self.inputs.get('lf_dark'))
        lf_live = store.get(self.inputs.get('lf_live'))
        mf = store.get(self.inputs.get('mf_live'))
        hf = store.get(self.inputs.get('hf'))
        drift = store.get(self.inputs.get('drift'))
        jumps = store.get(self.inputs.get('jumps'))

        results = {}

        results.update(lf_dark)
        results.update(lf_live)
        results.update(mf)
        results.update(hf)
        results.update(drift)
        results.update(jumps)

        self.logger.info("results: %s" % results)

        
