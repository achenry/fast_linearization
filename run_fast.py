import numpy as np
from weis.aeroelasticse.CaseGen_General import CaseGen_General
from weis.aeroelasticse.runFAST_pywrapper import runFAST_pywrapper_batch
from weis.aeroelasticse.LinearFAST import LinearFAST
import os
import multiprocessing as mp

ALL_DOFs = ['FlapDOF1', 'FlapDOF2', 'EdgeDOF', 'TeetDOF', 'DrTrDOF', 'GenDOF',
            'YawDOF', 'TwFADOF1', 'TwFADOF2', 'TwSSDOF1', 'TwSSDOF2',
            'PtfmSgDOF', 'PtfmSwDOF', 'PtfmHvDOF', 'PtfmRDOF', 'PtfmPDOF', 'PtfmYDOF']

ALL_COMPS = ['CompElast', 'CompInflow', 'CompAero', 'CompServo', 'CompHydro', 'CompSub', 'CompMooring', 'CompIce']

INIT_CONDITIONS = ['RotSpeed', 'BlPitch1', 'BlPitch2', 'BlPitch3', 'PtfmRoll', 'PtfmPitch', 'PtfmYaw', 'NacYaw']


def generate_linear_models(**kwargs):
    wind_indices = kwargs['wind_indices']

    case_inputs = {}

    for comp in ALL_COMPS:
        if comp == 'CompMooring':
            case_inputs[('Fst', comp)] = {'vals': [3 if comp in kwargs['comps'] else 0], 'group': 0}
        elif comp == 'CompAero':
            case_inputs[('Fst', comp)] = {'vals': [2 if comp in kwargs['comps'] else 0], 'group': 0}
        else:
            case_inputs[('Fst', comp)] = {'vals': [1 if comp in kwargs['comps'] else 0], 'group': 0}

    # elastodyn options
    for dof in ALL_DOFs:
        case_inputs[("ElastoDyn", dof)] = {'vals': ['True' if dof in kwargs['DOFs'] else 'False'
                                                    for case_idx in range(len(kwargs['wind_speeds']))], 'group': 1}

    case_list, case_name_list = CaseGen_General(case_inputs,
                                                dir_matrix=kwargs['fast_run_dir'],
                                                namebase=os.path.basename(kwargs['fast_model_dir']))

    fast_runner = LinearFAST(
        FAST_exe=kwargs['fast_exe_path'],
        FAST_lib=kwargs['fast_lib_path'],
        use_exe=True,
        FAST_InputFile=kwargs['fast_input_filename'],
        FAST_directory=kwargs['fast_model_dir'],
        FAST_runDirectory=kwargs['fast_run_dir'],
        cores=mp.cpu_count(),
        DT=kwargs['dt'],
        vrated='',
        wind_speeds=list(np.array(kwargs['wind_speeds'])[wind_indices]),
        TMax=kwargs['Tmax'],
        NLinTimes=36,
        DOFs=kwargs['DOFs'],
        wind_files=kwargs['wind_files'],
        case_list=case_list,
        case_name_list=case_name_list,
        TrimCase=3)

    fast_runner.cores = mp.cpu_count()
    fast_runner.FAST_exe = kwargs['fast_exe_path']
    fast_runner.FAST_lib = kwargs['fast_lib_path']
    fast_runner.FAST_InputFile = kwargs['fast_input_filename']
    fast_runner.FAST_directory = kwargs['fast_model_dir']
    fast_runner.FAST_runDirectory = kwargs['fast_run_dir']
    fast_runner.wind_speeds = list(np.array(kwargs['wind_speeds'])[wind_indices])
    fast_runner.TMax = kwargs['Tmax']
    fast_runner.NLinTimes = 36
    fast_runner.DOFs = kwargs['DOFs']
    fast_runner.rosco_infile = kwargs['rosco_discon_path']
    fast_runner.DT = kwargs['dt']
    fast_runner.use_exe = True
    fast_runner.rated_offset = 0.6
    fast_runner.v_rated = 11.4
    if 'discon_in' not in kwargs:
        kwargs['discon_in'] = {'VS_Rgn2K': None, 'PC_RefSpd': None, 'VS_RtTq': None}
    # fast_runner.fst_vt = {'DISCON_in': kwargs['discon_in']}

    fast_runner.gen_linear_cases(inputs={**kwargs['operating_points'], 'Vrated': 11.4, 'discon_in': kwargs['discon_in']})
    fast_runner.gen_linear_model()

    return fast_runner


def run_nonlin_simulations(**kwargs):
    # generate the matrix of cases
    case_inputs = {}

    wind_indices = kwargs['wind_indices']
    comps = kwargs['comps']

    # .fst options
    case_inputs[('Fst', 'TMax')] = {'vals': [kwargs['Tmax']], 'group': 0}
    case_inputs[('Fst', 'DT')] = {'vals': [kwargs['dt']], 'group': 0}
    case_inputs[('Fst', 'CalcSteady')] = {'vals': ['True'], 'group': 0}

    # servodyn options
    case_inputs[("ServoDyn", "DLL_FileName")] = {'vals': [kwargs['rosco_dll_path']], 'group': 0}
    case_inputs[("ServoDyn", "DLL_InFile")] = {'vals': [kwargs['rosco_discon_path']], 'group': 0}
    case_inputs[("ServoDyn", "PCMode")] = {'vals': [5], 'group': 0}
    case_inputs[("ServoDyn", "VSContrl")] = {'vals': [5], 'group': 0}
    case_inputs[("ServoDyn", "YCMode")] = {'vals': [0], 'group': 0}
    case_inputs[("ServoDyn", "AfCMode")] = {'vals': [0], 'group': 0}
    case_inputs[("ServoDyn", "CCMode")] = {'vals': [0], 'group': 0}
    case_inputs[("ServoDyn", "HSSBrMode")] = {'vals': [0], 'group': 0}

    if 'CompHydro' in comps:
        case_inputs[("HydroDyn", "WaveMod")] = {'vals': [0], 'group': 0}
        case_inputs[("HydroDyn", "RdtnMod")] = {'vals': [0], 'group': 0}
        case_inputs[("HydroDyn", "ExctnMod")] = {'vals': [0], 'group': 0}
        case_inputs[("HydroDyn", "PotMod")] = {'vals': [0], 'group': 0}

    for comp in ALL_COMPS:
        if comp == 'CompMooring':
            case_inputs[('Fst', comp)] = {'vals': [3 if comp in kwargs['comps'] else 0], 'group': 0}
        elif comp == 'CompAero':
            case_inputs[('Fst', comp)] = {'vals': [2 if comp in kwargs['comps'] else 0], 'group': 0}
        else:
            case_inputs[('Fst', comp)] = {'vals': [1 if comp in kwargs['comps'] else 0], 'group': 0}

    # # elastodyn options
    for dof in ALL_DOFs:
        case_inputs[("ElastoDyn", dof)] = {'vals': ['True' if dof in kwargs['DOFs'] else 'False'
                                                    for case_idx in range(len(wind_indices))],
        'group': 1}

    # single step wind inflow run
    if kwargs['wind_type'] == 'step':

        # set initial conditions
        case_inputs[("ElastoDyn", "RotSpeed")] = {'vals': [0],'group': 0}
        case_inputs[("ElastoDyn", "BlPitch1")] = {'vals': [0], 'group': 0}
        case_inputs[("ElastoDyn", "BlPitch2")] = {'vals': [0], 'group': 0}
        case_inputs[("ElastoDyn", "BlPitch3")] = {'vals': [0], 'group': 0}

        # set windtype to uniform and set input filename
        case_inputs[("InflowWind", "WindType")] = {'vals': [2], 'group': 0}
        case_inputs[("InflowWind", "Filename_Uni")] = {'vals': kwargs['wind_files'], 'group': 0}

        # case_inputs[("ServoDyn", "PCMode")] = {'vals': [5], 'group': 0}

    # multiple steady wind speeds
    elif kwargs['wind_type'] == 'steady':

        # set windtype to steady wind and set wind speeds
        case_inputs[("InflowWind", "WindType")] = {'vals': [1], 'group': 0}
        case_inputs[("InflowWind", "HWindSpeed")] = {'vals': [kwargs['wind_speeds'][i] for i in wind_indices], 'group': 1}

    # multiple uniform windfiles
    elif kwargs['wind_type'] == 'uniform':

        # set windtype to uniform wind and set wind files
        case_inputs[("InflowWind", "WindType")] = {'vals': [2], 'group': 0}
        case_inputs[("InflowWind", "Filename_Uni")] = {'vals': kwargs['wind_files'], 'group': 1}


    # set initial conditions
    if 'init_conditions' in kwargs:
        for key, value in  kwargs['init_conditions']:
            case_inputs[('ElastoDyn', key)] = {'vals': value[wind_indices], 'group': 1}


    # if kwargs['ol']:
    #     # if open-loop mode, supply open loop file with gentq and blpitch values
    #     case_inputs[("DISCON_in", "OL_Mode")] = {'vals': [1], 'group': 0}
    #     case_inputs[("DISCON_in", "OL_Filename")] = {'vals': kwargs['ol_filenames'], 'group': 1}
    #     case_inputs[("DISCON_in", "Ind_Breakpoint")] = {'vals': [1], 'group': 0}
    #     case_inputs[("DISCON_in", "Ind_BldPitch")] = {'vals': [2], 'group': 0}
    #     case_inputs[("DISCON_in", "Ind_GenTq")] = {'vals': [3], 'group': 0}
    #     case_inputs[("ServoDyn", "PCMode")] = {'vals': [0], 'group': 0}
    #
    # else:
    #     # else, for closed-loop mode, allow generator torque controller
    #     case_inputs[("DISCON_in", "OL_Mode")] = {'vals': [0], 'group': 0}

    case_list, case_name_list = CaseGen_General(case_inputs,
                                                dir_matrix=kwargs['fast_run_dir'],
                                                namebase="{model_name}_{wind_type}".format(
                                                    model_name=os.path.basename(kwargs['fast_input_filename']).split('.')[0],
                                                    wind_type=kwargs['wind_type']))

    fast_runner = runFAST_pywrapper_batch()

    fast_runner.FAST_exe = kwargs['fast_exe_path']
    fast_runner.FAST_lib = kwargs['fast_lib_path']
    fast_runner.FAST_InputFile = kwargs['fast_input_filename']
    fast_runner.FAST_directory = kwargs['fast_model_dir']
    fast_runner.FAST_runDirectory = kwargs['fast_run_dir']
    fast_runner.case_list = case_list
    fast_runner.case_name_list = case_name_list
    fast_runner.use_exe = True

    # fast_runner.run_multi()
    fast_runner.run_serial()
    return fast_runner
