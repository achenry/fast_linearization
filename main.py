# LINEARIZATION BUGS
# disturbances read from file MUST have same precision as u_ops, or else lin system simulation becomes unstable
# Do I need to set each blade pitch command individually as well as collective blade pitch command?
# y_ops are in deg, rpm, x_ops are in rad, rad/s, could also be the case for Nm for kNm ?
# setting PtfmSurge initial condition=> :Small angle assumption violated in SUBROUTINE SmllRotTrans() due to a large UR_bar input angles, SkewedWakeCorrection encountered a large value of chi,  BEM solution is being turned off due to low TSR.
# setting PtfmPitch => The BEM solution is being turned off due to low TSR., : SkewedWakeCorrection encountered a large value of chi

# LINEARIZATION TESTING NEXT STEPS
# Regenerate linear models without MoorDyn, SubDyn, HydroDyn and test?

# DEBUG NONLIN TODO
# test with Model_lin from command line with closed loop SKIPPED
# if work, openfast compilation is good
#   test with Model_lin from command line with rosco open loop ONGOING
#   if work, rosco open loop is good
#       test with Model_lin with rosco open loop from weis w/ no parameter changes ONGOING = MATCHES COMMAND LINE UP TO 29seconds TESTING TAB1
#       if work, weis is good
#           test with parameter changes one at a time
#           if work, param updates are good
#           else problem with param updates
#       else problem with weis
#   else problem with rosco open loop
# else problem is with openfast compilation

# ROSCO COMPILATION
# either use ifort, setting ifort alias with explicit library in bash_profile
# or use gfortran, setting SDKROOT env variable explicitly in bash_profile, and DON'T setvars.sh for ifort
# run setup.py install in WEIS to use ifort, or compile directly (cmake) with gfortran

# Questions
# What can we set DT to? => unchanged for Fst, 0.0025 for MoorDyn
# Which USFLOWT model can I use? => any
# Should T of .wnd files = T of simulation? What units of time are .wnd files in? Is DIST_DT written anywhere in model? => seconds, if final T of wind < Tmax of sim, final value holds
# What lb, ub, int do you recommend for step wind? => u0 - 0.5, u0 + 0.5, int=0.01?
# What amp, freq do you recommend for sine wind? => 0.5, freq?
# ED init conditions for blade pitch from outputs or inputs of linear model? => input (output should be identical, can't init output)
# ED where to find 'TeetDefl','TTDspFA', 'TTDspSS'  in linear model? => if disabled in ED then won't be in linear model, bending mode not the same as displacement
#  side-side to 0
# Do we run an OL simulation? => Turn control off in both or on in both.

# 1) Using  the LinearTurbineModel addControl and solve methods to simulate a linear system for perturbations in wind, pitch control or other inputs
# OL control => make a new controller class. use weis rather than ROSCO,
# aeroelasticsse/openmdao_openfast.py = main driver file, passes inputs to openfast and reads outputs from openfast

# 2) Running nonlinear OpenFAST simulations for user-given blade pitch, wind (or other input) time series using WEIS tools (or is this only possible in Simulink?)
# through ROSCO
# sim from ROSCO runs nonlin simulation

# 3) generating DLC wind cases
# modelling_options.yaml - DLC_driver, DLCGenerator
# ROSCO_testing.py (a bit outdated), ROSCO_testing_lite

# weis inputs: geometry.yaml (sim to fast input deck), analysis_options.yaml and modelling_options.yaml sets up weis
# openmdao - optimization framework

# TODO program ROSCO to compute steady-state dynamically and pause execution then
# TODO test sine sweep or superposition of sines

import control
import weis
import ROSCO_toolbox
import pyFAST
from ROSCO_toolbox.linear.linear_models import LinearTurbineModel
# from weis.control.LinearModel import LinearTurbineModel
# from ROSCO_toolbox.turbine import Turbine
from ROSCO_toolbox.control_interface import ControllerInterface as ROSCO_ControllerInterface

from weis.aeroelasticse.runFAST_pywrapper import runFAST_pywrapper_batch, libext
from weis.aeroelasticse.CaseGen_General import CaseGen_General
# from ROSCO_toolbox.ofTools.case_gen.CaseGen_IEC import CaseGen_IEC
# from ROSCO_toolbox.ofTools.fast_io.pyIECWind import pyIECWind_extreme, pyIECWind_turb
# from pyFAST.input_output.fast_linearization_file import FASTLinearizationFile
from pyFAST.input_output.fast_output_file import FASTOutputFile
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import os
import platform
import glob
import pickle

from scipy.interpolate import interp1d

from disturbances import generate_uniform_step_input, generate_uniform_sine_input, read_dist, generate_uniform_wind, plot_disturbances
import argparse

import multiprocessing as mp

from scipy.signal import lti, dlsim

np.seterr(all="raise")



def run_linear_simulation_wrapper(case, linturb_models, nonlin_output, dt, uni_tmax=0):
    if False:
        if case[("InflowWind", 'WindType')] == 2:
            wind = read_dist(case[("InflowWind", "Filename_Uni")])
        elif case[("InflowWind", 'WindType')] == 1:
            time = np.arange(0, TMAX, DT)
            wind = np.hstack(
                [time[:, np.newaxis], case[("InflowWind", 'HWindSpeed')] * np.ones(len(time))[:, np.newaxis]])
        disturbance['Time'] = wind[:, 0]
        disturbance['Wind'] = wind[:, 1]
        # H0.00-W0.0_2_DISCON.IN

        controller = ROSCO_ControllerInterface(lib_name=rosco_lib_path,
                                               param_filename=os.path.join('./results',
                                                                           f'{case_name}_DISCON.IN'))
        lin_outdata, lin_outlist, lin_plant = linturb_models.solve(disturbance=disturbance,
                                                                   controller=controller, open_loop=False)
        return {'outdata': pd.DataFrame(lin_outdata), 'lin_plant': lin_plant}
    else:
        dist = read_dist(case[("DISCON_in", "OL_Filename")])
        dist_time = dist[:, 0]
        # dist_blpitch = dist[:, 1]
        # dist_gentq = dist[:, 2]

        gentq_idx = linturb_models.DescCntrlInpt.index('ED Generator torque, Nm')
        # TODO set individual blade pitch commands too?
        blpitch_idx = [
            linturb_models.DescCntrlInpt.index(
                'ED Extended input: collective blade-pitch command, rad'),
            # linturb_models.DescCntrlInpt.index(
            #     'ED Blade 1 pitch command, rad'),
            # linturb_models.DescCntrlInpt.index(
            #     'ED Blade 2 pitch command, rad'),
            # linturb_models.DescCntrlInpt.index(
            #     'ED Blade 3 pitch command, rad')
        ]

        output = {channel: [] for channel in linturb_models.DescOutput}
        output['Time'] = []
        out_dt = nonlin_output['Time_[s]'].iloc[1] - nonlin_output['Time_[s]'].iloc[0]
        dt = 0.001
        time_ts = np.arange(uni_tmax, np.max(dist_time), dt)#np.array(nonlin_output['Time_[s]'])

        dist_u = interp1d(dist_time, dist[:, 1:], axis=0)(time_ts)

        wind_idx = list(linturb_models.u_h).index(case[('InflowWind', 'HWindSpeed')])

        u_ts = np.vstack([linturb_models.u_ops[:, wind_idx] for i in range(len(time_ts))])
        for i in blpitch_idx:
            u_ts[:, i] = dist_u[:, 0]
        u_ts[:, gentq_idx] = dist_u[:, 1]

        x, y = run_linear_simulations(linturb_models.x_ops[:, wind_idx],
                                      linturb_models.u_ops[:, wind_idx],
                                      linturb_models.y_ops[:, wind_idx],
                                      linturb_models.A_ops[:, :, wind_idx],
                                      linturb_models.B_ops[:, :, wind_idx],
                                      linturb_models.C_ops[:, :, wind_idx],
                                      linturb_models.D_ops[:, :, wind_idx],
                                      time_ts, u_ts, dt, out_dt=out_dt)

        for channel_idx, channel in enumerate(linturb_models.DescOutput):
            # output[channel] = np.concatenate([output[channel], y[:, channel_idx]])
            output[channel] = output[channel] + list(y[:, channel_idx])

        output['Time'] = np.concatenate([output['Time'], time_ts])

        return pd.DataFrame(output)


# def plot_case_grp(case_grp, n_disturbance_types, n_quantities, dist_cases, nonlin_outputs, lin_outputs):
def plot_case_grp(case_grp_df, nonlin_output_cols):
    n_disturbance_types = len(pd.unique(case_grp_df.index['DistType']))
    quantities = pd.unique(next(iter(case_grp_df)).index['Quantity'])
    n_quantities = len(quantities)

    # for each case e.g. each amplitude / frequency etc
    # for dist_case_idx, dist_case in dist_cases.iterrows():
    # for each DistSignal, DistType, HWindSpeed group
    for grp_idx, grp_df in case_grp_df:

        # initialize subplots
        fig, ax = plt.subplots(nrows=n_disturbance_types + n_quantities, ncols=1, sharex=True)
        ax[0].set_ylabel('BlPitch', rotation=0)
        ax[1].set_ylabel('GenTq', rotation=0)

        for q in range(n_quantities):
            ax[2 + q].set_ylabel(nonlin_output_cols[1 + q], rotation=0)

        ax[-1].set_xlabel('Time')

        dist_by = grp_df.groupby('OL_Filename')

        for dist_idx, dist_df in dist_by:
            # read in the disturbance signal
            dist_filename = dist_idx['OL_Filename']
            dist_data = read_dist(dist_filename)
            time = dist_data[:, 0]
            blpitch = dist_data[:, 1]
            gentq = dist_data[:, 2]

            # plot blade pitch signal
            ax[0].plot(time, blpitch, label=dist_idx)

            # plot generator torque signal
            ax[1].plot(time, gentq, label=dist_idx)

        # for each output quantity
        for q, quant in enumerate(pd.unique(dist_df['Quantity'])):
            # plot linear and nonlinear time series for each case of each quantity
            ax[2 + q].plot(time, nonlin_outputs[c].loc[nonlin_outputs.case_name == case.case_name, quant])
    plt.show()


def inv_mbc(dq_quants):
    T_inv = lambda theta: \
        np.vstack([
            [1, np.cos(theta), np.sin(theta)],
            [1, np.cos(theta + (2 * np.pi / 3)), np.sin(theta + (2 * np.pi / 3))],
            [1, np.cos(theta + (4 * np.pi / 3)), np.sin(theta + (4 * np.pi / 3))]
        ])
    return np.matmul(T_inv, dq_quants)


def setup_nonlinear_simulations(fast_model_dir, fast_install_dir, fast_dir, fast_input_file, Tmax, dt,
                                x_init, u_init, y_init,
                                path2dll, path2discon, dist_files, wind_files, debug=False):
    # instantiate run wrapper
    fast_exe = os.path.join(fast_install_dir, 'bin/openfast')
    print(fast_exe)
    fast_lib_dir = os.path.join(fast_install_dir, 'lib/libopenfastlib' + libext)
    print(fast_lib_dir)

    fast_runner = runFAST_pywrapper_batch(
        FAST_exe=fast_exe,
        FAST_lib=fast_lib_dir,
        FAST_runDirectory='./results',
        FAST_directory=fast_dir,
        FAST_InputFile=fast_input_file,
        use_exe=True)

    # User settings

    # settings passed to OpenFAST
    case_inputs = {}
    case_inputs[('Fst', 'TMax')] = {'vals': [Tmax], 'group': 0}
    # case_inputs[('Fst', 'DT')] = {'vals': [dt], 'group': 0}

    # NOTE: all Comp values must be turned on in files initially if they are to be used at all
    # NOTE: Replace Gonzalís with Gonzalis in Aerodyn.data
    # NOTE: Replace rosco.dat in ServoDyn.dat with rosco_ol.dat
    # NOTE: Replace DISCON with "DISCON" in ServoDyn.dat
    # NOTE: Change WindType to 2
    # NOTE: Change CpCqCt to "CpCqCt" in rosco.dat
    # NOTE: Change FATAL to "FATAL" in Model.fst

    # case_inputs[('Fst', 'CompElast')] = {'vals': [1], 'group': 0}
    # case_inputs[('Fst', 'CompInflow')] = {'vals': [1], 'group': 0}
    # case_inputs[('Fst', 'CompAero')] = {'vals': [2], 'group': 0}
    # case_inputs[('Fst', 'CompServo')] = {'vals': [1], 'group': 0}
    # case_inputs[('Fst', 'CompHydro')] = {'vals': [1], 'group': 0}
    # case_inputs[('Fst', 'CompSub')] = {'vals': [1], 'group': 0}
    # case_inputs[('Fst', 'CompMooring')] = {'vals': [3], 'group': 0}
    # case_inputs[('Fst', 'CompIce')] = {'vals': [0], 'group': 0}
    # case_inputs[('Fst', 'Linearize')] = {'vals': [False], 'group': 0}

    case_inputs[("AeroDyn15", "AFAeroMod")] = {'vals': [1], 'group': 0} # steady aerodynamics
    case_inputs[("AeroDyn15", "WakeMod")] = {'vals': [1], 'group': 0} # TODO check w/ David, also AeroDyn/SkewMod

    case_inputs[("HydroDyn", "WaveMod")] = {'vals': [0], 'group': 0} # turn off wave input
    case_inputs[("HydroDyn", "ExctnMod")] = {'vals': [0], 'group': 0} # turn off wave excitation
    case_inputs[("HydroDyn", "RdtnMod")] = {'vals': [0], 'group': 0} # turn off radiation memory effect, how platform influences waves in water

    # todo do we need some kind of control other than nonliner to keep from becoming unstable?
    # case_inputs[("ServoDyn", "DLL_FileName")] = {'vals': [path2dll], 'group': 0}
    case_inputs[("ServoDyn", "DLL_InFile")] = {'vals': [path2discon], 'group': 0}
    # case_inputs[("ServoDyn", "PCMode")] = {'vals': [5], 'group': 0}
    # case_inputs[("ServoDyn", "VSContrl")] = {'vals': [5], 'group': 0}
    # case_inputs[("ServoDyn", "YCMode")] = {'vals': [0], 'group': 0}
    # case_inputs[("ServoDyn", "AfCmode")] = {'vals': [0], 'group': 0}
    # case_inputs[("ServoDyn", "CCmode")] = {'vals': [0], 'group': 0}

    # blade pitch and generator torque disturbances
    case_inputs[("DISCON_in", "OL_Mode")] = {'vals': [1], 'group': 0}
    case_inputs[("DISCON_in", "Ind_Breakpoint")] = {'vals': [1], 'group': 0}
    case_inputs[("DISCON_in", "Ind_BldPitch")] = {'vals': [2], 'group': 0}
    case_inputs[("DISCON_in", "Ind_GenTq")] = {'vals': [3], 'group': 0}
    case_inputs[("DISCON_in", "PerfFileName")] = {'vals': [os.path.join(fast_dir, 'CpCtCq.dat')], 'group': 0}
    # TESTING TAB 3 => RUNNING NORMALLY
    case_inputs[("DISCON_in", "OL_Filename")] = {'vals': dist_files['filename'], 'group': 1}
    case_inputs[("DISCON_in", "DistType")] = {'vals': dist_files['dist_type'], 'group': 1}
    case_inputs[("DISCON_in", "DistSignal")] = {'vals': dist_files['dist_signal'], 'group': 1}

    n_dist_cases = len(dist_files['filename'])
    n_wind_cases = len(wind_files['speeds'])
    wind_repmat = lambda input: np.concatenate([
        [f for i in range(int(n_dist_cases / n_wind_cases))] for f in input
    ])
    case_inputs[("InflowWind", "WindType")] = {'vals': [1],
                                               'group': 0}  # 1=steady wind conditions, 2=uniform from file # TODO 1
    # case_inputs[("InflowWind", "Filename_Uni")] = {'vals': wind_repmat(wind_files['filenames']), 'group': 1}
    case_inputs[("InflowWind", "HWindSpeed")] = {'vals': wind_repmat(wind_files['speeds']), 'group': 1}

    # ElastoDyn initial conditions for each wind speed
    # for param in ['OoPDefl', 'IPDefl', 'BlPitch1', 'BlPitch2', 'BlPitch3', 'TeetDefl', 'Azimuth', 'RotSpeed', 'NacYaw',
    #               'TTDspFA', 'TTDspSS', 'PtfmSurge', 'PtfmSway', 'PtfmHeave', 'PtfmRoll', 'PtfmPitch', 'PtfmYaw']:

    # set initial conditions, converting rad to deg and rad/s to rpm
    # TESTING TAB 2 => STUCK FOR ALL USED,
    # => SD_CalcOutput: Angles in GetSmllRotAngs() are larger
    #  than 0.4 radians.
    # => The BEM solution is being turned off due to low TSR.
    # => Warning: SkewedWakeCorrection encountered a large value of chi (93.568 deg), so the yaw
    #  correction will be limited.
    # =>  Small angle assumption violated in SUBROUTINE SmllRotTrans() due to a large UR_bar input angles.
    case_inputs[('ElastoDyn', 'BlPitch1')] = {'vals': wind_repmat(
        u_init['ED Extended input: collective blade-pitch command, rad'] * RAD_TO_DEG
    ), 'group': 1}
    case_inputs[('ElastoDyn', 'BlPitch2')] = {'vals': wind_repmat(
        u_init['ED Extended input: collective blade-pitch command, rad'] * RAD_TO_DEG
    ), 'group': 1}
    case_inputs[('ElastoDyn', 'BlPitch3')] = {'vals': wind_repmat(
        u_init['ED Extended input: collective blade-pitch command, rad'] * RAD_TO_DEG
    ), 'group': 1}
    # # TODO QUESTION under MBC transform, do we just set Azimuth to 0?
    case_inputs[('ElastoDyn', 'Azimuth')] = {'vals': wind_repmat(
        np.zeros(n_wind_cases),
    ), 'group': 1}

    case_inputs[('ElastoDyn', 'RotSpeed')] = {'vals': wind_repmat(
        x_init[
            'ED First time derivative of Variable speed generator DOF (internal DOF index = DOF_GeAz), rad/s'] * RADS_TO_RPM
    ), 'group': 1}
    try:
        case_inputs[('ElastoDyn', 'NacYaw')] = {'vals': wind_repmat(
            x_init['ED Nacelle yaw DOF (internal DOF index = DOF_Yaw), rad'] * RAD_TO_DEG
        ), 'group': 1}
    except KeyError as e:
        print(e)
    # # # TODO does this correspond to correct state?
    case_inputs[('ElastoDyn', 'TTDspFA')] = {'vals': wind_repmat(
        x_init['ED 1st tower fore-aft bending mode DOF (internal DOF index = DOF_TFA1), m']
    ), 'group': 1}
    try:
        case_inputs[('ElastoDyn', 'TTDspSS')] = {'vals': wind_repmat(
            x_init['ED 1st tower side-to-side bending mode DOF (internal DOF index = DOF_TSS1), m']
        ), 'group': 1}
    except KeyError as e:
        print(e)

    # TESTING TAB1
    #      no run
    # y_init['ED PtfmSurge, (m)'] = 9.3255
    # case_inputs[('ElastoDyn', 'PtfmSurge')] = {'vals': wind_repmat(
    #     x_init['ED Platform horizontal surge translation DOF (internal DOF index = DOF_Sg), m']), 'group': 1}
    try:
        case_inputs[('ElastoDyn', 'PtfmSway')] = {'vals': wind_repmat(
            x_init['ED Platform horizontal sway translation DOF (internal DOF index = DOF_Sw), m']), 'group': 1}
    except KeyError as e:
        print(e)
    case_inputs[('ElastoDyn', 'PtfmHeave')] = {'vals': wind_repmat(
        x_init['ED Platform vertical heave translation DOF (internal DOF index = DOF_Hv), m']), 'group': 1}
    # y_init['ED PtfmRoll, (rad)']
    try:
        case_inputs[('ElastoDyn', 'PtfmRoll')] = {'vals': wind_repmat(
            x_init['ED Platform roll tilt rotation DOF (internal DOF index = DOF_R), rad'] * RAD_TO_DEG
        ), 'group': 1}
    except KeyError as e:
        print(e)
    # y_init['ED PtfmPitch, (rad)'] = 4.249 deg

    #      no run
    # case_inputs[('ElastoDyn', 'PtfmPitch')] = {'vals': wind_repmat(
    #     x_init['ED Platform pitch tilt rotation DOF (internal DOF index = DOF_P), rad'] * RAD_TO_DEG
    # ), 'group': 1}
    # y_init['ED PtfmYaw, (rad)']
    try:
        case_inputs[('ElastoDyn', 'PtfmYaw')] = {'vals': wind_repmat(
            x_init['ED Platform yaw rotation DOF (internal DOF index = DOF_Y), rad'] * RAD_TO_DEG
        ), 'group': 1}
    except KeyError as e:
        print(e)

    # generate the matrix of cases
    case_list, case_name_list = CaseGen_General(case_inputs,
                                                dir_matrix=fast_runner.FAST_runDirectory,
                                                namebase=os.path.basename(fast_dir))

    fast_runner.case_list = case_list
    fast_runner.case_name_list = case_name_list

    if PARALLEL:
        fast_runner.setup_multi()
    else:
        fast_runner.setup_serial()

    return fast_runner


def run_linear_simulations(x0, u0, y0, A, B, C, D, time_ts, u_ts, dt, out_dt=0.1):
    tout = []
    yout = []
    xout = []
    # # abs_x = x0
    delta_x = np.zeros_like(x0)
    # x_prev = x0

    # BUG: very unstable to small numerical errors
    system = lti(A, B, C, D).to_discrete(dt)

    # BUG precision issues arise when using absolute values rather than deviations
    if False:
        tout, yout, xout = dlsim(system, u_ts, time_ts, x0)
    else:
        A_d = system.A
        B_d = system.B
        # TODO tidy up
        for k in range(len(time_ts)):

            delta_u = u_ts[k] - u0

            try:
                delta_x = np.matmul(A_d, delta_x) + np.matmul(B_d, delta_u)
                delta_y = np.matmul(C, delta_x) + np.matmul(D, delta_u)
            except RuntimeWarning as e:
                print(e)
            except FloatingPointError as e:
                print(e)

            if k % (int(out_dt // dt)) == 0:
                tout.append(k * dt)
                xout.append(x0 + delta_x)
                yout.append(y0 + delta_y)

            if k % (int(1 // dt) + 1) == 0:
                print(f'Time = {k * dt} s, RotSpeed = {yout[-1][9]} rad/s, GenTq = {yout[-1][3]} N-m, BldPitch = {yout[-1][5]} rad')

    return tout, np.vstack(xout), np.vstack(yout)


def load_lin_models(linfile_root, parallel, wind_speeds=None):
    # if linfile_dirs is None:
    #     linfile_dirs = sorted([os.path.join(linfile_root, d) for d in os.listdir(linfile_root) if
    #                            os.path.isdir(os.path.join(linfile_root, d))])
    # if len(linfile_dirs) > 1:
    # os.listdir(linfile_root)
    lin_filenames = glob.glob(os.path.join(linfile_root, '**', '*.lin'), recursive=True)
    # elif len(linfile_dirs) == 1 and os.path.exists(os.path.join(linfile_root, linfile_dirs[0])):
    #     linfile_dirs = os.path.join(linfile_root, linfile_dirs[0])
    #     lin_filenames = glob.glob(os.path.join(linfile_dirs, '*.lin'))
    # else:
    #     linfile_dirs = linfile_root
    #     lin_filenames = glob.glob(os.path.join(linfile_root, '*.lin'))

    lin_files = [os.path.splitext(os.path.basename(file))[0].split('.') for file in
                 lin_filenames]  # [wind_speed, lin_number]
    lin_fileroots = sorted(set([f[0] for f in lin_files]))
    lin_numbers = set([f[1] for f in lin_files])

    linturb_models = LinearTurbineModel(linfile_root, lin_fileroots,
                                        nlin=max([int(num) for num in lin_numbers]), load_parallel=parallel)

    if wind_speeds is not None:
        wind_speed_idx = [list(linturb_models.u_h).index(u) for u in wind_speeds if
                          list(linturb_models.u_h).index(u) > -1]
        linturb_models.u_h = linturb_models.u_h[wind_speed_idx]
        linturb_models.A_ops = linturb_models.A_ops[:, :, wind_speed_idx]
        linturb_models.B_ops = linturb_models.B_ops[:, :, wind_speed_idx]
        linturb_models.C_ops = linturb_models.C_ops[:, :, wind_speed_idx]
        linturb_models.D_ops = linturb_models.D_ops[:, :, wind_speed_idx]
        linturb_models.x_ops = linturb_models.x_ops[:, wind_speed_idx]
        linturb_models.u_ops = linturb_models.u_ops[:, wind_speed_idx]
        linturb_models.y_ops = linturb_models.y_ops[:, wind_speed_idx]

    linturb_models.save(os.path.join(linfile_root, 'linturb_models'))

    return linturb_models


if __name__ == '__main__':

    fig_path = './results/figs/'
    if not os.path.exists(fig_path):
        os.mkdir(fig_path)

    parser = argparse.ArgumentParser()
    arg_flags = [('-p-', '--parallel'), ('-gd', '--gen_dist'), ('-tlm', '--trans_lin_models'),
                 ('-sns', '--setup_nonlin_sims'),
                 ('-rns', '--run_nonlin_sims'), ('-rls', '--run_lin_sims'), ('-pp', '--post_process'),
                 ('-d', '--debug')]
    for arg_flag in arg_flags:
        parser.add_argument(arg_flag[0], arg_flag[1], action="store_true")
    args = parser.parse_args()

    ## SETUP FLAGS, CONSTS, DIRECTORIES
    PARALLEL = args.parallel
    GENERATE_DISTURBANCES = args.gen_dist
    TRANSFORM_LINEAR_MODELS = args.trans_lin_models
    SETUP_NONLINEAR_SIMULATIONS = args.setup_nonlin_sims
    RUN_NONLINEAR_SIMULATIONS = args.run_nonlin_sims
    RUN_LINEAR_SIMULATIONS = args.run_lin_sims
    POST_PROCESS = args.post_process
    DEBUG = args.debug

    DT = 0.005
    TMAX = 5 if DEBUG else 1200
    WIND_SPEEDS = [20] #[12, 14, 20]

    ONE_DEGREE = 2 * np.pi / 360
    THOU_NM = 1000
    RPM_TO_RADS = 2 * np.pi / 60
    RADS_TO_RPM = 1 / RPM_TO_RADS
    DEG_TO_RAD = np.pi / 180
    RAD_TO_DEG = 1 / DEG_TO_RAD

    DIST_DT = 0.1
    UNI_TMAX = 1 if DEBUG else 600

    NONLIN_OUTPUT_COLS = ['Time_[s]', 'RotSpeed_[rpm]']  # , 'PtfmHeave_[m]', 'PtfmPitch_[deg]', 'PtfmSurge_[m]']
    LIN_OUTPUT_COLS = ['ED RotSpeed, (rad/s)']  # , 'ED PtfmHeave, (m)', 'ED PtfmPitch, (rad)', 'ED PtfmSurge, (m)']

    weis_dir = os.path.abspath(weis.__file__)
    rosco_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(ROSCO_toolbox.__file__)), 'ROSCO'))

    mactype = platform.system().lower()
    print(mactype)
    if mactype in ["linux", "linux2"]:
        # CHECK up-to-date models have been transferred to RC
        # scp -r /Users/aoifework/Documents/models/OpenFAST-Dev3.0_Models_Dist_v4 aohe7145@login.rc.colorado.edu:/projects/aohe7145/models/OpenFAST-Dev3.0_Models_Dist_v4
        # CHECK up-to-date weis toolbox
        # git add . && git commit -m "updats" && git push
        # git pull
        fast_install_dir = '/projects/aohe7145/src/openfast_usflowt/install'
        fast_model_dir = '/projects/aohe7145/models/OpenFAST-Dev3.0_Models_Dist_v4'
        model_name = 'H0.00-W0.0'
        # model_name = 'Hs1.13-WT11.4_C3'
        fast_input_file = 'Model.fst'
        linfile_root = '/projects/aohe7145/linearization/linearized_models'

        print(fast_install_dir)
        print(fast_model_dir)
        fast_dir = os.path.abspath(os.path.join(fast_model_dir, model_name))
        discon_path = os.path.abspath('../Model/Servo/rosco_ol.dat')

    elif mactype == "darwin":
        if False:  # requires weis-dev-env
            fast_install_dir = '/Users/aoifework/Documents/toolboxes/WEIS/local'
            model_name = 'IEA-15-240-RWT-UMaineSemi'
            fast_model_dir = f'/Users/aoifework/Documents/toolboxes/WEIS/ROSCO/Test_Cases/{model_name}'
            linfile_root = f'{fast_model_dir}/linearizations/'
            fast_input_file = 'IEA-15-240-RWT-UMaineSemi_ol.fst'
            fast_dir = os.path.abspath(fast_model_dir)
            discon_path = os.path.join(fast_model_dir, 'ServoData', 'DISCON-UMaineSemi_ol')
        else:  # requires weis-usflowt-env
            fast_install_dir = '/Users/aoifework/Documents/src/openfast_usflowt/install'
            fast_model_dir = '/Users/aoifework/Documents/models/OpenFAST-Dev3.0_Models_Dist_v4'
            model_name = 'H0.00-W0.0'
            linfile_root = '/Users/aoifework/Documents/linearization/linearized_models'
            fast_dir = os.path.abspath(os.path.join(fast_model_dir, model_name))
            discon_path = "/Users/Documents/models/OpenFAST-Dev3.0_Models_Dist_v4/Model/Servo/rosco_ol.dat"
            fast_input_file = 'Model.fst'

            # fast_model_dir = '/Users/aoifework/Documents/models/Model_lin'
            # model_name = 'Hs1.13-WT11.4_C3'
            # fast_input_file = 'Model_U_14.fst'

    rosco_lib_path = os.path.join(rosco_dir, 'install/lib/libdiscon' + libext)
    print(rosco_lib_path)

    ## SETUP LINEAR SIMULATIONS
    if TRANSFORM_LINEAR_MODELS:
        linturb_models = load_lin_models(linfile_root, parallel=PARALLEL,
                                         wind_speeds=WIND_SPEEDS)
    else:
        linturb_models = LinearTurbineModel.load(os.path.join(linfile_root, 'linturb_models'))

    # check stability and controllability
    if False:
        for u_idx in range(len(linturb_models.u_h)):
            A = linturb_models.A_ops[:, :, u_idx]
            B = linturb_models.B_ops[:, :, u_idx]
            ctrb_mat = control.ctrb(A, B)
            is_ctrb = np.linalg.matrix_rank(ctrb_mat) == ctrb_mat.shape[0]
            eig_vals = np.linalg.eigvals(A)
            is_stable = np.all(np.real(eig_vals) < 0)
            print(f'stability = {is_stable} for u={linturb_models.u_h[u_idx]}')
            print(f'controllability = {is_ctrb} for u={linturb_models.u_h[u_idx]}')

    ## FETCH EQUILIBRIUM POINTS FROM LINEAR MODELS

    # for lin_q_idx, lin_q in enumerate(linturb_models.DescOutput):
    #     if 'rad/s' in lin_q:
    #         linturb_models.y_ops[lin_q_idx] = linturb_models.y_ops[lin_q_idx] * RPM_TO_RADS
    #     elif 'rad' in lin_q:
    #         linturb_models.y_ops[lin_q_idx] = linturb_models.y_ops[lin_q_idx] * DEG_TO_RAD

    x_desc = linturb_models.DescStates
    x_init = {k: v for (k, v) in zip(x_desc, linturb_models.x_ops)}
    u_desc = linturb_models.DescCntrlInpt
    u_init = {k: v for (k, v) in zip(u_desc, linturb_models.u_ops)}
    y_desc = linturb_models.DescOutput
    y_init = {k: v for (k, v) in zip(y_desc, linturb_models.y_ops)}

    uni_wind_files = {'speeds': [], 'filenames': []}
    for u_idx, u in enumerate(linturb_models.u_h):
        filedir = os.path.join(fast_model_dir, 'Model', 'Wind')
        if not os.path.exists(filedir):
            os.makedirs(filedir)
        filename = os.path.join(filedir, f'uni-wind-{u}.dat')
        uni_wind_files['filenames'].append(filename)
        uni_wind_files['speeds'].append(u)
        generate_uniform_wind(filename, DIST_DT, TMAX, [u])

    ## GENERATE DISTURBANCES
    # for each input in lin_data from IfW of ED
    # generate a step and sinusoidal series

    # Ask Nikhar re units of control signals from ROSCO => rad and Nm
    dist_inc = 0.001
    blpitch_inc = dist_inc * u_init['ED Extended input: collective blade-pitch command, rad']
    gentq_inc = dist_inc * u_init['ED Generator torque, Nm']
    gentq_dist_mult = np.array([0, -1, 1])
    # blpitch_dist_mult = np.array([])
    blpitch_dist_mult = np.array([0, -1, 1])
    dist_freq = [0.1]
    # dist_freq = np.array([0.1, 0.2, 0.4, 1.4, 2])

    disturbances = {
        'pitch': {  # in rad
            'uni_val': u_init['ED Extended input: collective blade-pitch command, rad'],
            'step_vals': np.array([blpitch_inc[u_idx] * blpitch_dist_mult for u_idx in range(len(linturb_models.u_h))]),
            'amp_vals': [[blpitch_inc[u_idx]] for u_idx in range(len(linturb_models.u_h))],
            # 'freq_vals': np.array([[] for u in linturb_models.u_h])
            'freq_vals': np.array([dist_freq for u in linturb_models.u_h])
        },
        'torque': {  # in Nm
            'uni_val': u_init['ED Generator torque, Nm'],
            'step_vals': np.array([gentq_inc[u_idx] * gentq_dist_mult for u_idx in range(len(linturb_models.u_h))]),
            'amp_vals': [[gentq_inc[u_idx]] for u_idx in range(len(linturb_models.u_h))],
            'freq_vals': np.array([dist_freq for u in linturb_models.u_h])  # np.array([0.1, 0.2, 0.4, 1.4, 2])
        }
    }

    dist_files = {'filename': [], 'wind_speed': [], 'dist_type': [], 'dist_signal': [],
                  'step_val': [], 'amp_val': [], 'freq_val': []}

    if GENERATE_DISTURBANCES:
        # for each wind
        for wind_idx, wind in enumerate(linturb_models.u_h):
            uni_vals = [disturbances[dist_type]['uni_val'][wind_idx] for dist_type in disturbances]

            # for each disturbance
            for dist_idx, dist_type in enumerate(disturbances):

                # set step_val, amp_val to zero for all other disturbances
                step_vals = [0 for dist in disturbances]
                amp_vals = [0 for dist in disturbances]
                freq_vals = [0 for dist in disturbances]

                # generate uniform (disturbance free) cases
                # uni_input_filename = os.path.join(fast_model_dir,
                #                                   f'Model/Servo/uni-{dist_type}.dat')
                # uni_input = generate_uniform_input(filename=uni_input_filename, dt=DIST_DT, Tmax=TMAX,
                #                                    uniform_val=uni_val)

                # vary step_val, amp_val, freq_val for this disturbance

                # generate step disturbance cases
                for step_idx, step_val in enumerate(disturbances[dist_type]['step_vals'][wind_idx]):
                    varied_step_vals = list(step_vals)
                    varied_step_vals[dist_idx] = step_val
                    step_dir = 'pos' if step_val >= 0 else 'neg'

                    filedir = os.path.join(fast_model_dir, 'Model', 'Servo')
                    if not os.path.exists(filedir):
                        os.makedirs(filedir)
                    filename = os.path.join(filedir, f'uni-step-{dist_type}-{step_dir}_{step_idx}-wind_{wind_idx}.dat')
                    uni_step_input = generate_uniform_step_input(filename, DIST_DT, UNI_TMAX, TMAX,
                                                                 uniform_vals=uni_vals,
                                                                 step_vals=varied_step_vals)

                    dist_files['filename'].append(filename)
                    dist_files['wind_speed'].append(wind)
                    dist_files['dist_type'].append(dist_type)
                    dist_files['dist_signal'].append('step_up' if step_val >= 0 else 'step_down')
                    dist_files['step_val'].append(abs(step_val))
                    dist_files['amp_val'].append(None)
                    dist_files['freq_val'].append(None)

                fixed_amp_vals = list(amp_vals)
                amp_idx = 0
                fixed_amp_vals[dist_idx] = disturbances[dist_type]['amp_vals'][wind_idx][amp_idx]
                amp_val = fixed_amp_vals[dist_idx]

                for freq_idx, freq_val in enumerate(disturbances[dist_type]['freq_vals'][wind_idx]):

                    varied_freq_vals = list(freq_vals)
                    varied_freq_vals[dist_idx] = freq_val

                    filedir = os.path.join(fast_model_dir, 'Model', 'Servo')
                    if not os.path.exists(filedir):
                        os.makedirs(filedir)
                    filename = os.path.join(filedir,
                                            f'uni-sine-{dist_type}-amp_{amp_idx}-freq_{freq_idx}-wind_{wind_idx}.dat')
                    uni_sine_input = generate_uniform_sine_input(filename, DIST_DT, UNI_TMAX, TMAX,
                                                                 uniform_vals=uni_vals,
                                                                 amp_vals=fixed_amp_vals,
                                                                 freq_vals=varied_freq_vals,
                                                                 phase_vals=[0, 0])

                    dist_files['filename'].append(filename)
                    dist_files['wind_speed'].append(wind)
                    dist_files['dist_type'].append(dist_type)
                    dist_files['dist_signal'].append('sine_freq')
                    dist_files['step_val'].append(None)
                    dist_files['amp_val'].append(amp_val)
                    dist_files['freq_val'].append(freq_val)

        pickle.dump(dist_files, open(os.path.join('./results', 'dist_files'), "wb"))
    else:
        dist_files = pickle.load(open(os.path.join('./results', 'dist_files'), 'rb'))


    PLOT_DISTURBANCES = False
    if PLOT_DISTURBANCES:
        plot_disturbances(disturbances, dist_files, fig_path)

    ## SETUP NONLINEAR SIMULATIONS
    if SETUP_NONLINEAR_SIMULATIONS:
        fast_runner = setup_nonlinear_simulations(fast_model_dir=fast_model_dir,
                                                  fast_install_dir=fast_install_dir,
                                                  fast_dir=fast_dir,
                                                  fast_input_file=fast_input_file,
                                                  Tmax=TMAX, dt=DT,
                                                  x_init=x_init, u_init=u_init, y_init=y_init,
                                                  path2dll=rosco_lib_path,
                                                  path2discon=discon_path,
                                                  dist_files=dist_files,
                                                  wind_files=uni_wind_files,
                                                  debug=DEBUG)
        pickle.dump(fast_runner, open(os.path.join('./results', 'fast_runner'), "wb"))
    else:
        fast_runner = pickle.load(open(os.path.join('./results', 'fast_runner'), 'rb'))

    ## RUN NONLINEAR SIMULATIONS
    if RUN_NONLINEAR_SIMULATIONS:
        if PARALLEL:
            fast_runner.run_multi()
        else:
            fast_runner.run_serial()
            # nonlin_outputs[0].loc[:, NONLIN_OUTPUT_COLS] # in-plane

        if PARALLEL:
            pool = mp.Pool(processes=mp.cpu_count())
            nonlin_outputs = pool.starmap(read_outfile, zip(fast_runner.case_name_list,
                                                            [fast_runner for c in fast_runner.case_name_list]))
        else:
            nonlin_outputs = []
            for c_name in fast_runner.case_name_list:
                nonlin_outputs.append(read_outfile(c_name, fast_runner))

        pickle.dump(nonlin_outputs, open(os.path.join(fast_runner.FAST_runDirectory, 'nonlin_outputs'), "wb"))

    else:
        nonlin_outputs = pickle.load(open(os.path.join(fast_runner.FAST_runDirectory, 'nonlin_outputs'), 'rb'))

    ## RUN LINEAR SIMULATIONS
    lin_outlist = linturb_models.getOutList()
    if RUN_LINEAR_SIMULATIONS:
        disturbance = {}

        if PARALLEL:
            pool = mp.Pool(processes=mp.cpu_count())
            lin_outputs = pool.starmap(run_linear_simulation_wrapper, zip(fast_runner.case_list,
                                                                          [linturb_models
                                                                           for c in fast_runner.case_list],
                                                                          nonlin_outputs,
                                                                          [DT for c in fast_runner.case_list],
                                                                          [UNI_TMAX for c in fast_runner.case_list]))
        else:
            lin_outputs = []
            for case_idx, case in enumerate(fast_runner.case_list):
                lin_outputs.append(run_linear_simulation_wrapper(case, linturb_models, nonlin_outputs[case_idx], DT, UNI_TMAX))

        pickle.dump(lin_outputs, open(os.path.join(fast_runner.FAST_runDirectory, 'lin_outputs'), "wb"))
    else:
        lin_outputs = pickle.load(open(os.path.join(fast_runner.FAST_runDirectory, 'lin_outputs'), 'rb'))

    # COMPARE NONLINEAR TO LINEAR RESULTS
    # steady-state error - range over the last 30sec < 10% of meter
    if POST_PROCESS:

        # fetch nonlinear results

        # ss_error = {'CaseName': [], 'DistFilename': [], 'DistType': [], 'DistSignal': [],
        #             'Quantity': [], 'HWindSpeed': [], 'SteadyStateError': []}

        # generate dataframes of linear vs nonlinear time-series

        case_results = list(fast_runner.case_list)
        case_result_names = list(fast_runner.case_name_list)
        ts_dfs = {}

        # for each case and each output quantity, load the time-series of the linear and nonlinear simulations into dataframes
        for c, case in enumerate(case_results):
            case_name = case_result_names[c]
            time_ts = lin_outputs[c]['Time']
            # for each quantity
            NEW_NONLIN_OUTPUT_COLS = list(NONLIN_OUTPUT_COLS)
            for nonlin_q, lin_q in zip(NONLIN_OUTPUT_COLS[1:], LIN_OUTPUT_COLS):
                wind_speed = case[('InflowWind', 'HWindSpeed')]

                # get the nonlinear time series

                # check units and convert rpm to rad/s
                if 'rpm' in nonlin_q:
                    new_nonlin_q = nonlin_q.replace('rpm', 'rad/s')
                    NEW_NONLIN_OUTPUT_COLS[NONLIN_OUTPUT_COLS.index(nonlin_q)] = new_nonlin_q
                    nonlin_outputs[c][nonlin_q] = nonlin_outputs[c][nonlin_q] * RPM_TO_RADS
                    nonlin_outputs[c].rename(columns={nonlin_q: new_nonlin_q}, inplace=True)
                    nonlin_q = new_nonlin_q
                elif 'deg' in nonlin_q:
                    new_nonlin_q = nonlin_q.replace('deg', 'rad')
                    NEW_NONLIN_OUTPUT_COLS[NONLIN_OUTPUT_COLS.index(nonlin_q)] = new_nonlin_q
                    nonlin_outputs[c].rename(columns={nonlin_q: new_nonlin_q}, inplace=True)
                    nonlin_outputs[c][new_nonlin_q] = nonlin_outputs[c][new_nonlin_q] * DEG_TO_RAD
                    nonlin_q = new_nonlin_q

                nonlin_ts = nonlin_outputs[c][nonlin_q]

                # get the linear time series
                # lin_ts = lin_outputs[c].loc[lin_outputs[c]['IfW Wind1VelX, (m/s)'] == wind_speed][lin_q]
                lin_ts = lin_outputs[c][lin_q]

                # TODO transform rotating quantities via inverse MBC

                # set quantity label
                # case_results[c]['Quantity'] = nonlin_q
                case_results[c][f'TimeTS'] = time_ts
                case_results[c][f'NonlinearTS_{nonlin_q}'] = nonlin_ts
                case_results[c][f'LinearTS_{nonlin_q}'] = lin_ts
                case_results[c][f'AbsDiffTS_{nonlin_q}'] = np.abs(nonlin_ts - lin_ts)

                # set quantity steady state error
                # case_results['SteadyStateError'][c] = \
                #     np.abs(nonlin_quant.loc[time == ss_time] - lin_quant.loc[time == ss_time])

                # find steady-state time
                # ss_time = nonlin_outputs[c].loc[
                #     ((nonlin_quant.rolling(30).max() - nonlin_quant.rolling(30).min()).abs() < 232),
                #     NONLIN_OUTPUT_COLS[0]].iloc[0]

        NONLIN_OUTPUT_COLS = list(NEW_NONLIN_OUTPUT_COLS)

        # split into groups for each wind speed and iterate through
        case_results_df = pd.DataFrame(case_results)
        case_results_df['case_name'] = case_result_names
        case_results_df['case_idx'] = np.arange(len(case_result_names))
        # new_cols = ['HWindSpeed', 'DistType', 'DistSignal', 'OL_Filename']
        for col in case_results_df.columns:
            if type(col) is tuple:
                new_col = col[1]
                case_results_df.rename(columns={col: new_col}, inplace=True)

        # measure the absolute time-varying error
        #  by_case = case_results_df.groupby(['HWindSpeed', 'OL_Filename'])
        #  case, frame = next(iter(by_case))
        # (by_case['NonlinearTS'].values() /
        # - by_case['LinearTS'].values()).abs()

        # by_dist = case_results_df.groupby(['HWindSpeed', 'DistType', 'DistSignal'])
        # for case_grp, case_grp_df in by_dist:
        #     plot_case_grp(case_grp_df=case_grp_df)

        dist_types = pd.unique(pd.unique(case_results_df['DistType']))
        n_disturbance_types = len(disturbances.keys())  # or disturbances.keys()
        quantities = NONLIN_OUTPUT_COLS[1:]
        n_quantities = len(quantities)
        wind_speeds = pd.unique(case_results_df['HWindSpeed'])
        for wind_speed in wind_speeds:
            wind_speed_cases = case_results_df.loc[case_results_df['HWindSpeed'] == wind_speed]

            for dist_type in dist_types:
                dist_type_cases = wind_speed_cases.loc[wind_speed_cases['DistType'] == dist_type]
                dist_signals = pd.unique(dist_type_cases['DistSignal'])

                for dist_signal in dist_signals:
                    # range of different disturbances for given type (blade pitch or torque) and signal (step-up, step-down, sine-freq)
                    dist_cases = dist_type_cases.loc[dist_type_cases['DistSignal'] == dist_signal]

                    # plot the BlPitch, GenTq and output signals
                    # plot_case_grp(dist_cases)

                    # for each case e.g. each amplitude / frequency etc
                    for dist_case_idx, dist_case in dist_cases.iterrows():

                        # initialize subplots
                        fig, ax = plt.subplots(nrows=n_disturbance_types + n_quantities, ncols=1, sharex=True)
                        ax[0].set_title(f'Wind Speed = {wind_speed} m/s')
                        ax[0].set_ylabel('BlPitch', rotation=0)
                        # if n_disturbance_types == 2:
                        ax[1].set_ylabel('GenTq', rotation=0)

                        for q in range(n_quantities):
                            ax[n_disturbance_types + q].set_ylabel(quantities[q], rotation=0)

                        ax[-1].set_xlabel('Time')

                        # read in the disturbance signal
                        dist_filename = dist_case['OL_Filename']
                        dist_data = read_dist(dist_filename)
                        dist_time = dist_data[:, 0]
                        blpitch = dist_data[:, 1]
                        gentq = dist_data[:, 2]

                        # plot blade pitch signal
                        ax[0].plot(dist_time, blpitch, label=f"Case {dist_case['case_idx']}")

                        # plot generator torque signal
                        # if n_disturbance_types == 2:
                        ax[1].plot(dist_time, gentq, label=f"Case {dist_case['case_idx']}")

                        # for each output quantity
                        for q, quant in enumerate(quantities):
                            nonlin_ts = dist_case[f'NonlinearTS_{quant}']
                            lin_ts = dist_case[f'LinearTS_{quant}']
                            error_ts = dist_case[f'AbsDiffTS_{quant}']

                            # plot linear and nonlinear time series for each case of each quantity
                            # if this works can remove time-series data from case list
                            # nonlin_ts = nonlin_outputs[c].loc[nonlin_outputs[c].case_name == dist_case.case_name, quant]
                            # lin_ts = lin_outputs[c].loc[nonlin_outputs.case_name == case.case_name, quant]
                            #     error_ts = np.abs(nonlin_ts - lin_ts)
                            time_ts = dist_case[f'TimeTS']
                            # ax[n_disturbance_types + q].plot(time_ts, nonlin_ts, label='Nonlinear')
                            ax[n_disturbance_types + q].plot(time_ts, lin_ts, label='Linear')
                            # ax[n_disturbance_types + q].plot(time_ts, error_ts, label='Error')
                            ax[n_disturbance_types + q].legend(loc='center right')

                        plt.show()
                        plt.savefig(os.path.join(fig_path, f'ts-case_{dist_case.case_idx}.png'))

                        # # add ss error as annotation
                        # ax[2 + q].annotate(f'ss error = {case.SteadyStateError}')

                        # for val in ['step_val', 'amp_val', 'freq_val']:

                    # plot the range of the disturbance signal over all cases on the top 2 subplots

                # plot the steady-state error for each quantity on a separate plot over the same range of disturbance signals

            # for c, case in case_results_df.iterrows():

            # dist_type = case[("DISCON_in", "DistType")]
            # dist_signal = case[("DISCON_in", "DistSignal")]
