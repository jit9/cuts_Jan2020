from todloop import TODLoop

from cuts import CutSources
loop = TODLoop()

loop.add_tod_list('./data/test.txt')

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

loop.run(0,1)
    

