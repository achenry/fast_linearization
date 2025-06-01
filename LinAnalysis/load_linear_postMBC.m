function [MBC, matData, FAST_linData, sys] = load_linear_postMBC(outdir, uref, outputs, rm_hydro)
% Loads outputs from openfast linearizations (10/6/2020)

% Inputs:
%     outdir:       str, directory of .lin files
%     OutfileBase:  str, base name of .lin files
%     nlin:         int, number of .lin files per operating point
%     outputs:      cell, desired outputs from linear models
%     rm_hydro:     int, (1/0 -> T/F) - remove hydrodynamic states


%% Define FileNames
FAST_linData_filename = [outdir '/' num2str(uref) '_FAST_lin.mat'];
matData_filename = [outdir '/' num2str(uref) '_matData.mat'];
MBC_filename = [outdir '/' num2str(uref) '_mbc_data.mat'];


%% Multiblade Coordinate Transform
% Since MBC has already been performed, just load the data

FAST_linData = load(FAST_linData_filename, 'FAST_linData').FAST_linData;
matData = load(matData_filename, 'matData').matData;
MBC = load(MBC_filename, 'mbc_data').mbc_data;

%% Remove inputs, states, and outputs

% Find input indices to keep
wind_inp     = find(contains(FAST_linData(1).u_desc, 'IfW Extended input: horizontal wind speed'));
cpit_inp     = find(contains(FAST_linData(1).u_desc, 'ED Extended input: collective blade-pitch command'));
torq_inp     = find(contains(FAST_linData(1).u_desc, 'ED Generator torque'));
inp_inds = [wind_inp, cpit_inp, torq_inp]
inpnames = FAST_linData(1).u_desc(inp_inds);

% Find state indices to remove
if rm_hydro == 1
    hydro_states = find(contains(MBC.DescStates,'HD'));
else
    hydro_states = [];
end
az_state     = find(contains(MBC.DescStates,'ED Variable speed generator DOF'));
twr_state    = []; find(contains(MBC.DescStates,'tower'));
heave_state  = []; find(contains(MBC.DescStates,'heave'));
surge_state  = []; find(contains(MBC.DescStates,'surge'));

rmstates = [hydro_states; az_state; twr_state; heave_state; surge_state];
statenames = MBC.DescStates;

% findoutput states to keep
output_inds = find(contains(FAST_linData(1).y_desc,outputs));
outnames =  FAST_linData(1).y_desc(output_inds);

% Rename for simplification
A_MBC = MBC.AvgA;
B_MBC = MBC.AvgB; 
C_MBC = MBC.AvgC;
D_MBC = MBC.AvgD;

% Conversions
rpm2rads = [find(contains(FAST_linData(1).y_desc,'rpm'))];
deg2rad  = [find(contains(FAST_linData(1).y_desc,'deg'))];
C_MBC(rpm2rads,:) = C_MBC(rpm2rads,:) * pi/30;
C_MBC(deg2rad,:)  = C_MBC(deg2rad,:) * pi/180;

% Keep inputs
B_MBC = B_MBC(:,inp_inds);
D_MBC = D_MBC(:,inp_inds);

% Remove states
A_MBC(rmstates,:) = []; A_MBC(:,rmstates) = [];
B_MBC(rmstates,:) = []; 
statenames(rmstates) = [];
C_MBC(:,rmstates) = [];
% B_MBC(:,rminputs) = [];
% D_MBC(:,rminputs) = [];

% Keep outputs
C_MBC = C_MBC(output_inds,:);
D_MBC = D_MBC(output_inds,:);

% Save matrices
MBC.AvgA = A_MBC;
MBC.AvgB = B_MBC;
MBC.AvgC = C_MBC;
MBC.AvgD = D_MBC;
MBC.DescStates  = statenames;
MBC.DescInps    = inpnames;
MBC.DescOutputs = outnames;

% Save state smpace
sys = ss(A_MBC, B_MBC, C_MBC, D_MBC,...
         'StateName',statenames, ...
         'InputName',inpnames, ...
         'OutputName', outnames);


% Save operating points
MBC.ops.wind   = FAST_linData(1).u_op{wind_inp}
MBC.ops.pitch  = FAST_linData(1).u_op{cpit_inp}
MBC.ops.torque = FAST_linData(1).u_op{torq_inp}

