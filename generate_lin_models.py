# PROBLEM could be with openfast_usflowt OR rosco
# could I update

# FIXES
# AeroDyn.dat gonzales
# update rosco.dat with missing lines
# Model.fst FATAL => "FATAL"
# ServoDyn DISCON -> "DISCON"
# MoorDyn WaterKin.dat -> "WaterKin.dat"

# QUESTIONS
# Should I turn off Hydro too?
# How should I set the controller (ServoDyn) for the nonlinear runs?
# How should other params be set for nonlinear runs?
# How to find rated wind speed?
# Which initial conditions need to be set to generate the linear models?
# How do I set fixed platform?
# try with Hs1.13WT11.4
# Cpcqct file in discon, DT, value of initial RotSpeed, using uni wind file as opposed to step, increased ServoDyn values

# Turning off PC_Node, setting blpitch init conds, try high eg 20/30 for 12, 40/50 for 20, should be stable
# running steady wind


# linearized for 20m/s
# plot outputs from failed runs
# try setting init conditions
# compare exported fast input files to Model_lin

# TESTING
# RUN1: Step 12-20
# RUN@: parallele 12, 14, 20, setting init conditions

# Package original C1 files and rewritten files and send to David

# from pyFAST.linearization.linearization import writeLinearizationFiles

from weis.aeroelasticse.runFAST_pywrapper import mactype, libext
import multiprocessing as mp
from disturbances import generate_uniform_input, generate_wind_files
from read_outputs import read_outfile
import pickle
from collections import defaultdict
import os
import shutil
import matplotlib.pyplot as plt
import numpy as np
from run_fast import generate_linear_models, run_nonlin_simulations
from scipy.io import loadmat


def main():
    # fast_model_version = 'OpenFAST-Dev3.2_Models_Baseline'
    # fast_model_version = 'NREL-5MW-cases'
    fast_model_version = ''

    if mactype == "darwin":
        home_dir = '/Users/aoifework/Documents'
        # model_root_dir = '/Users/aoifework/Documents/dev/WEIS/ROSCO/Test_Cases'
        # model_root_dir = '/Users/aoifework/Documents/Research/transition_mpc/models/'
        model_root_dir = os.path.join(home_dir, 'Research/ipc_tuning/SOAR-25-V2f_IF')

        # openfast_root_dir = '/Users/aoifework/Documents/usflowt_src/openfast/'
        openfast_install_dir = '/Users/aoifework/Documents/toolboxes/openfast/install'
        # openfast_install_dir = os.path.join(home_dir, 'usflowt_src/openfast/install')
        # fast_models_dir = os.path.join(model_root_dir, fast_model_version)
        fast_models_dir = os.path.join(model_root_dir, fast_model_version)

        wind_dir = os.path.join(fast_models_dir, 'WindData')

    elif mactype in ["linux", "linux2"]:
        root_dir = '/projects/aohe7145/'

        fast_models_dir = os.path.join(root_dir, 'src')

    rosco_dll_path = os.path.join(openfast_install_dir, 'lib', 'libdiscon{ext}'.format(ext=libext))
    fast_exe_path = os.path.join(openfast_install_dir, 'bin', 'openfast')
    fast_lib_path = os.path.join(openfast_install_dir, 'lib', 'libopenfastlib{ext}'.format(ext=libext))
    fast_model_dir = fast_models_dir

    # fast_input_filename = 'NREL-5MW_aboveR.fst'
    fast_input_filename = 'weis_job_00.fst'

    rosco_discon_path = os.path.join(fast_model_dir, 'DISCON.IN')
    if not os.path.exists(rosco_discon_path):
        rosco_discon_path = None
        print('rosco_discon_path does not exist')
    discon_in = {
        'VS_Rgn2K': 1.5480e+08,
        # Generator torque constant in Region 2 (HSS side), [Nm/(rad/s)^2], approx 2.18575 from NREL 5W, 1.5480e+08 from k_opt_actual
        'PC_RefSpd': 0.5506,
        # Desired (reference) HSS speed for pitch controller, [rad/s], identical to VS_RefSpd usually, approx 122.90967 from NREL 5W, 0.5506 from Parameters.Turbine.wg_rated
        'VS_RtTq': 45403.4244e+03,
        # Rated torque, [Nm], approx 43093.51876 from NREL 5MW, 45403.4244e+03 from Parameters.Turbine.T_rated
        'WE_GearboxRatio': 1
    }

    # define run parameters
    WIND_TYPE = 'steady'
    # WIND_TYPE = 'step'
    OL = False
    wind_speeds = list(np.arange(0.5, 24.5, 0.5))  # [22.5]
    wind_indices = np.arange(len(wind_speeds))
    Tmax = 250
    dt = 0.0125

    linmodel_type = 'excGenDOF_incSecOrd'
    if 'incSecOrd' in linmodel_type:
        DOFs = ['FlapDOF1', 'FlapDOF2', 'EdgeDOF', 'TwFADOF1', 'TwFADOF2', 'TwSSDOF1', 'TwSSDOF2'] # 'GenDOF'
    else:
        DOFs = ['FlapDOF1', 'EdgeDOF', 'TwFADOF1', 'TwSSDOF1']

    comps = ['CompElast', 'CompInflow', 'CompAero', 'CompServo']  # , 'CompSub', 'CompMooring', 'CompHydro']

    run_type = '{wind_type}_wind-{control_mode}'.format(wind_type=WIND_TYPE, control_mode="OL" if OL else "CL")

    lin_model_version = '{version}_compall'.format(version=fast_model_version)

    fast_nonlin_run_dir = os.path.join(fast_model_dir,
                                       'linearization/{run_type}/nonlinear_results/{linmodel_type}'.format(
                                           run_type=run_type, linmodel_type=linmodel_type))
    if not os.path.exists(fast_nonlin_run_dir):
        os.makedirs(fast_nonlin_run_dir)

    fast_lin_run_dir = os.path.join(fast_model_dir, 'linearization/{run_type}/lin_results/{linmodel_type}'.format(
        run_type=run_type, linmodel_type=linmodel_type))
    if not os.path.exists(fast_lin_run_dir):
        os.makedirs(fast_lin_run_dir)

    fast_linfile_dir = os.path.join(fast_model_dir, 'linearization/{run_type}/linfiles/{linmodel_type}'.format(
        run_type=run_type, linmodel_type=linmodel_type))
    if not os.path.exists(fast_linfile_dir):
        os.makedirs(fast_linfile_dir)

    # DOFs = ['FlapDOF1',  'GenDOF']
    #          # 'PtfmSgDOF', 'PtfmSwDOF', 'PtfmHvDOF', 'PtfmRDOF', 'PtfmPDOF', 'PtfmYDOF']



    if WIND_TYPE in ['step', 'uniform']:
        wind_files = generate_wind_files(Tmax=Tmax, dt=dt, fast_models_dir=fast_models_dir,
                                         wind_speeds=wind_speeds[wind_indices],
                                         wind_type=WIND_TYPE, t_int=50, wind_dir=wind_dir, zero_start=False)
        generate_wind_files(Tmax=Tmax, dt=dt, fast_models_dir=fast_models_dir, wind_speeds=wind_speeds[wind_indices],
                            wind_type=WIND_TYPE, t_int=50, wind_dir=wind_dir, zero_start=False)
    else:
        wind_files = None

    # if using Open Loop control mode, generate time series input for blade pitch and generator torque
    if OL:

        # uni_vals for blpitch and gentq for wind_speeds
        uni_vals = {'GenTq': np.array([2.072E+05] * 3)[wind_indices],
                    'BlPitch': np.array([0.06421, 0.1577, 0.3073])[wind_indices]}

        ol_filenames = [os.path.join(fast_models_dir, 'Model', 'Servo', 'ol_inputs_{int(u)}.dat'.format(u=wind_speed))
                        for wind_speed in wind_speeds]
        for f_idx, f in enumerate(ol_filenames):
            generate_uniform_input(f, dt=dt, Tmax=Tmax,
                                   uniform_vals=[uni_vals['BlPitch'][f_idx], uni_vals['GenTq'][f_idx]])
    else:
        ol_filenames = None

    RUN_NONLIN = False

    if RUN_NONLIN:
        nonlin_fast_runner = run_nonlin_simulations(
            fast_models_dir=fast_models_dir, fast_model_dir=fast_model_dir, fast_exe_path=fast_exe_path,
            fast_lib_path=fast_lib_path, fast_run_dir=fast_nonlin_run_dir, fast_input_filename=fast_input_filename,
            wind_speeds=wind_speeds, rosco_dll_path=rosco_dll_path, rosco_discon_path=rosco_discon_path, Tmax=Tmax,
            dt=dt, wind_indices=wind_indices,
            wind_files=wind_files, DOFs=DOFs, comps=comps, wind_type=WIND_TYPE, ol=OL, ol_filenames=ol_filenames
        )

        pickle.dump(nonlin_fast_runner, open(os.path.join(fast_nonlin_run_dir, 'nonlin_fast_runner'), "wb"))

    elif os.path.exists(os.path.join(fast_nonlin_run_dir, 'nonlin_fast_runner')):
        nonlin_fast_runner = pickle.load(open(os.path.join(fast_nonlin_run_dir, 'nonlin_fast_runner'), 'rb'))

    READ_NONLIN = False
    if READ_NONLIN:
        # fetch steady-state values for initial conditions
        pool = mp.Pool(processes=mp.cpu_count())
        runDir_list = fast_nonlin_run_dir  #
        cname_list = ["{model_name}_{wind_type}_{c}".format(
            model_name=os.path.basename(fast_input_filename).split('.')[0],
            wind_type=WIND_TYPE, c=c) for c in range(len(wind_speeds))]
        nonlin_outputs = pool.starmap(read_outfile, zip(cname_list,
                                                        [nonlin_fast_runner.FAST_runDirectory for c in cname_list]))
        pickle.dump(nonlin_outputs, open(os.path.join(fast_nonlin_run_dir, 'nonlin_outputs'), "wb"))
    elif os.path.exists(os.path.join(fast_nonlin_run_dir, 'nonlin_outputs')):
        nonlin_outputs = pickle.load(open(os.path.join(fast_nonlin_run_dir, 'nonlin_outputs'), 'rb'))

    # fetch ss operating points
    op_point_key = {
        'OoPDefl': 'TipDxc1_[m]',
        'IPDefl': 'TipDyc1_[m]',
        'TeetDefl': None,
        'BlPitch1': 'BldPitch1_[deg]',
        # 'BlPitch2': 'BldPitch2_[deg]',
        # 'BlPitch3': 'BldPitch3_[deg]',
        'RotSpeed': 'RotSpeed_[rpm]',
        # 'NacYaw': 'NacYaw_[deg]',
        'TTDspFA': None,  # 'YawBrTDxt',
        'TTDspSS': None,  # 'YawBrTDyt',
        'GenTq': 'GenTq_[kN-m]'
        # 'PtfmSurge': 'PtfmSurge_[m]',
        # 'PtfmSway': 'PtfmSway_[m]',
        # 'PtfmHeave': 'PtfmHeave_[m]',
        # 'PtfmRoll': 'PtfmRoll_[deg]',
        # 'PtfmPitch': 'PtfmPitch_[deg]',
        # 'PtfmYaw': 'PtfmYaw_[deg]'
    }
    outputs = ['RotSpeed_[rpm]', 'BldPitch1_[deg]', 'RootMyb1_[kN-m]', 'GenTq_[kN-m]']
    # plot initial condition quantities
    # outputs = [quant for quant in op_point_key.values() if quant is not None]
    if 'nonlin_outputs' in locals():
        nonlin_outputs_gen = (op for i, op in enumerate(nonlin_outputs) if i in wind_indices)

        fig, ax = plt.subplots(len(outputs), 1)
        for case_idx, case in enumerate(nonlin_outputs_gen):
            for q_idx, q in enumerate(outputs):
                # print(case_idx)
                # print(case)
                ax[q_idx].plot(case['Time_[s]'].array, case[q].array,
                               label='{u} m/s'.format(u=np.array(wind_speeds)[case_idx]))
                ax[q_idx].set_title(q)
            # ax[-1].legend()
        if WIND_TYPE == 'steady' and (not OL):
            ax[0].set_title('Steady Wind Speed, PCMode off, DLL Gentq controller')
        elif OL:
            ax[0].set_title('Steady Wind Speed, PCMode off, OL Blpitch and Gentq commands')
        elif WIND_TYPE == 'step':
            ax[0].set_title('Step Wind Speeds, PCMode on, DLL Blpitch and Gentq commands')
        plt.show()

        nonlin_outputs_gen = (op for i, op in enumerate(nonlin_outputs) if i in wind_indices)
        operating_points = defaultdict(list)
        for case in nonlin_outputs_gen:
            for init_cond_key, q in op_point_key.items():
                if q is None or q not in case.columns:
                    continue
                operating_points[init_cond_key].append(case[q].rolling(window=1000).mean().iloc[-1])
    else:
        # read from mat file
        ss_vals = loadmat(os.path.join(model_root_dir, 'ss_vals.mat'))
        operating_points = defaultdict(list)
        for case_idx, case in enumerate(np.array(wind_speeds)[wind_indices]):
            wind_speed_idx = wind_speeds.index(case)
            for init_cond_key, q in op_point_key.items():
                if init_cond_key in ss_vals['ss_vals'][0].dtype.names:
                    operating_points[init_cond_key].append(ss_vals['ss_vals'][0][0][ss_vals['ss_vals'][0].dtype.names.index(init_cond_key)][0][wind_speed_idx])
                else:
                    operating_points[init_cond_key].append(0)


    # generate linear models @ 12, 14, 20 m/s
    GEN_LIN = True
    if GEN_LIN:
        lin_fast_runner = generate_linear_models(
            fast_models_dir=fast_models_dir, fast_model_dir=fast_model_dir, fast_exe_path=fast_exe_path,
            fast_lib_path=fast_lib_path, fast_run_dir=fast_lin_run_dir, fast_input_filename=fast_input_filename,
            wind_speeds=wind_speeds, wind_files=wind_files, wind_indices=wind_indices,
            rosco_dll_path=rosco_dll_path, rosco_discon_path=rosco_discon_path, Tmax=Tmax, dt=dt,
            DOFs=DOFs, operating_points=operating_points, fast_linfile_dir=fast_linfile_dir, comps=comps,
            discon_in=discon_in)

        # move linear models
        for file in [f for f in os.listdir(fast_lin_run_dir) if f.endswith('.lin')]:
            shutil.move(os.path.join(fast_lin_run_dir, file), os.path.join(fast_linfile_dir, file))

        # fetch steady-state values for initial conditions
        cname_list = [f'lin_{c:02}' for c in range(len(wind_indices))]  # lin_fast_runner.case_name_list
        runDir_list = [fast_lin_run_dir for c in range(len(wind_indices))]
        pool = mp.Pool(processes=mp.cpu_count())
        lin_outputs = pool.starmap(read_outfile, zip(cname_list, runDir_list))

        pickle.dump(lin_fast_runner, open(os.path.join(fast_lin_run_dir, 'lin_fast_runner'), "wb"))
        pickle.dump(lin_outputs, open(os.path.join(fast_lin_run_dir, 'lin_outputs'), "wb"))
    else:
        # lin_fast_runner = pickle.load(open(os.path.join(fast_lin_run_dir, 'lin_fast_runner'), 'rb'))
        lin_outputs = pickle.load(open(os.path.join(fast_lin_run_dir, 'lin_outputs'), 'rb'))

    # plot linearization trajectories
    fig, ax = plt.subplots(len(outputs), 1)
    # lin_outputs_gen = (op for i, op in enumerate(lin_outputs) if i in wind_indices)

    for case_idx, case in enumerate(lin_outputs):
        # if case['Time_[s]'].iloc[-1] < 240 / dt:
        #     print(case_idx)
        for q_idx, q in enumerate(outputs):
            ax[q_idx].plot(case['Time_[s]'].array, case[q].array, label='{u}'.format(u=case_idx))
            ax[q_idx].set_title(q)
    ax[0].legend()
    plt.show()


if __name__ == '__main__':
    main()
