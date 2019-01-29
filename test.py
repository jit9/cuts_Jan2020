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
    'output_key': 'scan_params',
}
loop.add_routine(AnalyzeScan(**scan_params))

thermal_params = {
    'input_key': 'tod',
    'output_key': 'thermal_results',
    'channel': None,
    'T_max': 0.10,
    'dT_max': 0.0015, 
}
loop.add_routine(AnalyzeTemperature(**thermal_params))

BASE_DIR = '/data/actpol/actpol_data_shared/ArrayData/2016/ar3/' 
gd_params = {
    'input_key': 'tod',
    'output_key': 'dets',
    'source': 'individual',
    'live': BASE_DIR + 'live_pa3_f90_s16_c10_v4.dict',
    'dark': BASE_DIR + 'dark.dict',
    'exclude': BASE_DIR + 'exclude_pa3_f90_s16_c10_v4.dict'
}
loop.add_routine(GetDetectors(**gd_params))

cal_params = {
    'input_key': 'tod',
    'dets_key': 'dets',
    'output_key': 'tod',
    'output_calData': 'calData',
    'flatfield': "/data/actpol/actpol_data_shared/FlatFields/2015/" + \
                 "ff_actpol3_2015_c9_w1_v2b_mix90-150_it9_actpol3_2015_c9_w2_photon_mix90-150_it7.dict",
    'config': [{
        "type": "depot_cal",
        "depot": "/data/actpol/depot",
        "tag": "pa3_s16_BS",
        "name": "biasstep"
    }, {
        "type": "constant",
        "value": 0.821018,
        "name": "DC bias factor"
    }],
    'forceNoResp': True,
    'calibrateTOD': True,
}
loop.add_routine(CalibrateTOD(**cal_params))

jump_params = {
    'input_key': 'tod',
    'output_key': 'crit_jumps',
    'dsStep': 4,
    'window': 1,
}
loop.add_routine(FindJumps(**jump_params))

fft_params = {
    'input_key': 'tod',
    'output_key': 'tod',
    'fft_data': 'fft_data'
}
loop.add_routine(FouriorTransform(**fft_params))

lf_dark_params = {
    'fft_data': 'fft_data',
    'dets': 'dets',
    'tod': 'tod',
    'output_key': 'lf_dark',
    'scan': 'scan_params',
    'presel': {
        'method': 'median',
        'minSel': 0,
        'initCorr': 0.9,
    },
    'useTaper': False,
    'cancelSync': False,
    'doubleMode': False,
    'freqRange' : {
        'fmin': 0.017,
        'fshift': 0.009,
        'band': 0.071,
        'Nwin': 1,
    },
}
loop.add_routine(AnalyzeDarkLF(**lf_dark_params))


# run pipeline
loop.run(100,101)
