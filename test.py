from todloop import TODLoop
from todloop.tod import TODLoader
from features import JesseFeatures

# initialize the pipeline
loop = TODLoop()

# specify the list of tods to go through
loop.add_tod_list("data/2016_ar3_train.txt")


################################
# add routines to the pipeline #
################################

# add a routine to load tod
loader_params = {
    'output_key': 'tod',
}
loop.add_routine(TODLoader(**loader_params))

params = {
    'inputs': {
        'tod': 'tod',
    },
    'outputs': {
        'results': 'jesse_features',
    }
}
loop.add_routine(JesseFeatures(**params))

# run the pipeline for the first pipeline
loop.run(0, 1)
