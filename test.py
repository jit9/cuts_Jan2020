from todloop import TODLoop
from todloop.tod import TODLoader

from cuts import *
loop = TODLoop()

loop.add_tod_list('/home/lmaurin/TODLists/2016_ar3_season_nohwp.txt')

loop.add_routine(TODLoader(output_key="tod", repair_pointing=False))

# I included all the parameters here but actually I don't need this
# much because all of the cuts have been generated before so I will
# only need to know the tag to retrieve it in principle.
source_params = {
    'input_key': 'tod',
    'output_key': 'tod',
    'tag_source': 'pa3_f90_s16_c10_v1_source',
    'source_list': '/home/lmaurin/TODLists/2016_ar3_season_nohwp.txt',
    'no_noise': True,
    'pointing_par': {
        'source': 'fp_file', \
        'filename': '/data/actpol/actpol_data_shared/RelativeOffsets/' + \
                    'template_ar3_s16_170131.txt'
    },
    'mask_shift_generator': {
    'source_list': '/data/actpol/actpol_data_shared/BrightSources/sources.txt', 
    'mask_params': { 
        'radius': (3./60)  #degrees
    },
    'mask_shift_generator': {\
        'source':'file',\
        'filename':'/data/actpol/actpol_data_shared/TODOffsets/' + \
                    'tod_offsets_2016_170131.txt',
        'columns': [0,3,4],
        'rescale_degrees': 1./60}
    },
    'depot': '/data/actpol/depot',
}
loop.add_routine(CutSources(**source_params))

planets_params = {
    'input_key': 'tod',
    'output_key': 'tod',
    'tag_planet': 'pa3_f90_s16_c10_v1_planet',
    'depot': '/data/actpol/depot',
}
loop.add_routine(CutPlanets(**planets_params))

sync_params = {
    'input_key': 'tod',
    'output_key': 'tod',
    'tag_sync': 'pa3_f90_s16_c10_v1',
    'remove_sync': False,
    'force_sync': False,
    'depot': '/data/actpol/depot',
}
loop.add_routine(RemoveSyncPickup(**sync_params))

partial_params = {
    'input_key': 'tod',
    'output_key': 'tod',
    'tag_partial': 'pa3_f90_s16_c10_v1_partial',
    'force_partial': False,
    'glitchp': { 'nSig': 10., 'tGlitch' : 0.007, 'minSeparation': 30, \
                 'maxGlitch': 50000, 'highPassFc': 6.0, 'buffer': 200 },
    'depot': '/data/actpol/depot',
}
loop.add_routine(CutPartial(**partial_params))

# hwp is not needed here, just added in for completeness
# hwp_params = {
#     'input_key': 'tod',
#     'output_key': 'tod',
#     'substract_hwp': False,
# }
# loop.add_routine(SubstractHWP(**hwp_params))

transform_params = {
    'remove_mean': False,
    'remove_median': True,
    'detrend': False,
    'remove_filter_gain': False
}
loop.add_routine(TransformTOD(**transform_params))

scan_params = {
    'input_key': 'tod',
    'output_key': 'chunk_params',
    'scan_param': {}
}
loop.add_routine(AnalyzeScan(**scan_params))

# run pipeline
loop.run(100,101)
    

