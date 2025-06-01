%% Generate CL Models
% 
% Load linearized models and a ROSCO gain schedule, generate closed-loop
% models from them. 

addpath(genpath('../../toolboxes/matlab-toolbox'));

%% ------- Linearized models -------
% -----------------------------------
%
% -- CASES -- :
% case_1: PtfmPDOF, GenDOF
% case_2: PtfmPDOF, GenDOF, TwFADOF1
% case_3: PtfmPDOF, GenDOF, TwFADOF1, PtfmSgDOF
% case_4: PtfmPDOF, GenDOF, TwFADOF1, PtfmSgDOF, PtfmHvDOF

% File setup
% outdir = 'iea15mw_cases/case_4';
outfiles_base = 'Model_U';
nlin = 12;              % Number of linearization points per wind speed

outdir = '../linearized_models';
postMBC_outdir = '../linearized_models';
uref_list = [12, 14, 20];%14:24;
outfiles_num = length(uref_list);      % Number of wind speeds the system was linearized at 

%% Get linearizattion outfile bases
for idx = 1:outfiles_num
     all_bases{idx} = [outfiles_base '_' num2str(uref_list(idx))];
end

%% Load all state space matrices
all_systems = {};
all_MBS     = {};
all_matdata = {};
all_linData = {};
all_sys     = {};

outputs = {'PtfmHeave', 'PtfmPitch', 'PtfmSurge', ...
            ...%'PtfmRoll', 'PtfmSway', 'PtfmYaw', ...
            'GenSpeed', 'GenPwr'};
        
for i = 1:length(uref_list)
    rm_hydro = 1;
    [MBC, matData, FAST_linData, sys] = load_linear(outdir, all_bases{i}, nlin, outputs, rm_hydro);
     %[MBC, matData, FAST_linData, sys] = load_linear_postMBC(postMBC_outdir, uref_list(i), outputs, rm_hydro);

    windspeeds(i) = MBC.WindSpeed;
    all_MBC{i} = MBC;
    all_matData{i} = matData;
    all_linData{i} = FAST_linData;
    all_sys{i} = sys;
end

%% Combine into 3d array for interpolation
A_set = zeros(length(windspeeds), size(all_MBC{i}.AvgA,1),  size(all_MBC{i}.AvgA,2));
B_set = zeros(length(windspeeds), size(all_MBC{i}.AvgB,1),  size(all_MBC{i}.AvgB,2));
C_set = zeros(length(windspeeds), size(all_MBC{i}.AvgC,1),  size(all_MBC{i}.AvgC,2));
D_set = zeros(length(windspeeds), size(all_MBC{i}.AvgD,1),  size(all_MBC{i}.AvgD,2));
pitch_op_set = zeros(1, length(windspeeds));
torque_op_set = zeros(1, length(windspeeds));

for i = 1:length(windspeeds)
    A_set(i,:,:) = all_MBC{i}.AvgA;
    B_set(i,:,:) = all_MBC{i}.AvgB;
    C_set(i,:,:) = all_MBC{i}.AvgC;
    D_set(i,:,:) = all_MBC{i}.AvgD;
    pitch_op_set(i) = all_MBC{i}.ops.pitch;
    torque_op_set(i) = all_MBC{i}.ops.torque;
end

%% Simulate OL linear models
t = 0:0.001:10;
for i = 1:1%length(windspeeds)
    sys = all_sys{i};
    % x0 = all_matData{i}.Avgxop;
    % step(sys);
    du = zeros(length(t), length(all_MBC{i}.DescInps));
    % pitch_dist = pitch_op_set(i) * 0.0001 * sin(t);
    dist = torque_op_set(i) * 0.001;% * sin(t);
%     input_string = 'collective blade-pitch';
    input_string = 'Generator torque';
%     % input_string = 'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s'
    inpt_idx = find(contains(all_MBC{i}.DescInps, input_string));
    du(:, inpt_idx) = dist;
    lsim(sys, du, t);
end

%% ------ Plant model for specific wind speed ------
% Interpolated linear models from openfast linearizations

% ws_des = 	16.5;        % Desired wind speed
% 
% A_op = interp1(windspeeds,A_set,ws_des);
% A_op = reshape(A_op, size(A_op,2), size(A_op,3));
% B_op = interp1(windspeeds,B_set,ws_des);
% B_op = reshape(B_op, size(B_op,2), size(B_op,3));
% C_op = interp1(windspeeds,C_set,ws_des);
% C_op = reshape(C_op, size(C_op,2), size(C_op,3));
% D_op = interp1(windspeeds,D_set,ws_des);
% D_op = reshape(D_op, size(D_op,2), size(D_op,3));
% 
% ss_op = ss(A_op, B_op, C_op, D_op, ...
%     'StateName', all_MBC{1}.DescStates, 'InputName', all_MBC{1}.DescInps, 'OutputName', all_MBC{1}.DescOutputs);
% 
% input_string = 'collective blade-pitch';
% % input_string = 'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s'
% uind = find(contains(all_MBC{1}.DescInps, input_string));
% 
% % output_string = 'ED PtfmPitch, (deg)';
% output_string = 'GenSpeed';
% yind = find(contains(all_MBC{1}.DescOutputs, output_string));
% 
% % Transfer Function
% Gs = tf(ss_op(yind,uind));


%% --------- Load Tuning Models ------------
% Load linear models from ROSCO tuning process (see example_11.py)
% rosco_model_file = './lin_models/IEA15MW_LinMod.dat'; % 
% rosco_model_file = './lin_models/IEA15MW_LinMod.dat'; % 
% [ModelData] = importdata(rosco_model_file,'\t',1);
% 
% for i = 1:length(ModelData.colheaders) % remove whitespace in column headers
%     ModelData.colheaders{i} = ModelData.colheaders{i}(ModelData.colheaders{i} ~= ' ');
% end
% 
% % Pull out tuning model values
% windspeed_tune  = ModelData.data(:,contains(ModelData.colheaders,'WindSpeed'));
% A_tune          = ModelData.data(:,contains(ModelData.colheaders,'A_om'));
% B_tau_tune      = ModelData.data(:,contains(ModelData.colheaders,'b_tau'));
% B_beta_tune     = ModelData.data(:,contains(ModelData.colheaders,'b_theta'));

%% -------- Dynamic rlocus, nyquist --------
% This generates dynamic root locus or nyquist plots. 
% For the Root Locus, S.fun should be @get_tf_cl or @get_tf
% For the Nyquist, S.fun should be @get_tf_ol
%
% The structure, S, is used to pass parameters to/from the function calls. 

% ws_rated = 11.4;
ws_rated = 14;

ws_0     = 16; ceil(ws_rated);
omega_0  = 0.3;
zeta_0   = 0.8;

omega_min = 0.05;
omega_max = 0.5;
zeta_min = 0.2;
zeta_max = 2.0;

% Setup figure
S.fh = figure('Units', 'normalized', 'Position', [0 0 1 1]);
% figure(S.fh, 'Units', 'normalized', 'Position', [0 0 1 1]);
S.ax = axes('Parent',S.fh,'units','normalized','position',[0.1 0.27 0.8 0.7]);

% Define function and plot type
S.fun    = @get_tf_cl;   % @get_tf_ol, @get_tf
S.pltfun = @rlocusplot;  % @nyquistplot

% Turbine and analysis space parameters
S.windspeeds = windspeeds;
S.all_MBC    = all_MBC;
S.ws_tune    = windspeed_tune;
S.A_tune     = A_tune;
S.B_tune     = B_beta_tune;
S.omega      = omega_0;
S.zeta       = zeta_0;
S.ws         = ws_0;

% Gain schedule parameters (from rosco_c4.in)
S.c4_oppoint = [0.076694  0.104632  0.126607  0.145560  0.162625  0.178367  0.193079  0.206934  0.220177  0.232974  0.245218  0.257054  0.268731  0.279796  0.290951  0.301449  0.312168  0.322206  0.332428  0.342278  0.351957  0.361895  0.371018  0.380343  0.389920  0.398608  0.407499];
S.c4_kp = [-0.028155  -0.024262  -0.021128  -0.018550  -0.016392  -0.014560  -0.012985  -0.011616  -0.010415  -0.009354  -0.008408  -0.007561  -0.006798  -0.006106  -0.005476  -0.004901  -0.004373  -0.003887  -0.003438  -0.003021  -0.002634  -0.002274  -0.001938  -0.001623  -0.001327  -0.001049  -0.000788];
S.c4_ki = [-0.005563  -0.004963  -0.004479  -0.004082  -0.003749  -0.003466  -0.003223  -0.003012  -0.002827  -0.002663  -0.002517  -0.002387  -0.002269  -0.002162  -0.002065  -0.001976  -0.001895  -0.001820  -0.001751  -0.001686  -0.001627  -0.001571  -0.001519  -0.001471  -0.001425  -0.001382  -0.001342];


% initial function eval and plot
S.Gs = S.fun(S);
S.h  = S.pltfun(S.Gs);
pfinfo = functions(S.pltfun);
if strcmp(pfinfo.function, 'nyquistplot')
    xlim = [-2, 3];
    ylim = [-6, 2];
elseif strcmp(pfinfo.function, 'rlocusplot')
%     xlim = [-20, 40];
%     ylim = [-50, 50];
    xlim = [-1, 1];
    ylim = [-10, 10];
end
S.pltopts.XLim = xlim;
S.pltopts.YLim = ylim;
options = getoptions(S.h);
options.XLim = S.pltopts.XLim;
options.YLim = S.pltopts.YLim;
update(S)

% Wind Speed Slider
S.slider1 = uicontrol('Parent',S.fh,'Style','slider',...
              'SliderStep',[0.1/(max(windspeeds)-ceil(ws_rated)), 0.1],...
              'Units','normalized','position',[0.1 0.18 0.8 0.02],...
              'value', S.ws, 'min', ceil(ws_rated), 'max', max(windspeeds),...
              'callback', {@SliderCB, 'ws'});
uicontrol('Parent',S.fh,'Style','text',...
              'Units','normalized','position',[0.1 0.15 0.8 0.03],... 
              'String','Wind Speed');
uicontrol('Parent',S.fh,'Style','text',...
              'Units','normalized','position',[0.08 0.14 0.04 0.03],... 
              'String',ceil(ws_rated)); 
uicontrol('Parent',S.fh,'Style','text',...
              'Units','normalized','position',[0.88 0.14 0.04 0.03],... 
              'String',max(windspeeds));

% If its not just the plant model, include control parameter sliders
if ~strcmp(char(S.fun),'get_tf')    
    S.slider2 = uicontrol('Parent',S.fh,'Style','slider',...
                  'SliderStep', [0.01/(omega_max-omega_min), 0.1],...
                  'Units','normalized','position',[0.1 0.11 0.8 0.02],...
                  'value', S.omega, 'min', omega_min, 'max', omega_max,...
                  'callback', {@SliderCB, 'omega'});
    uicontrol('Parent',S.fh,'Style','text',...
              'Units','normalized','position',[0.1 0.09 0.8 0.03],... 
              'String','omega_pc');
    uicontrol('Parent',S.fh,'Style','text',...
                  'Units','normalized','position',[0.08 0.09 0.04 0.02],... 
                  'String',omega_min);
    uicontrol('Parent',S.fh,'Style','text',...
                  'Units','normalized','position',[0.88 0.09 0.04 0.02],... 
                  'String',omega_max);   

    S.slider3 = uicontrol('Parent',S.fh,'Style','slider',...
                  'SliderStep', [0.01/(zeta_max-zeta_min), 0.1],...
                  'Units','normalized','position',[0.1 0.045 0.8 0.02],...
                  'value', S.zeta, 'min', zeta_min, 'max', zeta_max,...
                  'callback', {@SliderCB, 'zeta'});
    uicontrol('Parent',S.fh,'Style','text',...
              'Units','normalized','position',[0.1 0.02 0.8 0.03],... 
              'String','zeta_pc');
    uicontrol('Parent',S.fh,'Style','text',...
                  'Units','normalized','position',[0.08 0.025 0.04 0.02],... 
                  'String',zeta_min);
    uicontrol('Parent',S.fh,'Style','text',...
                  'Units','normalized','position',[0.88 0.025 0.04 0.02],... 
                  'String',zeta_max);   
end
guidata(S.fh, S);  % Store S structure in the figure


%% Sanity check
% S.omega = get(S.slider2, 'Value');
% S.ws = get(S.slider1, 'Value');
% 
% [Hs, Gs, Cs] = get_tf_cl(S);

%% Get plant transfer function
function [Gs, ss_op] = get_tf(S)
    %{ 
    Get plant transfer functions. Input is a structure that has the
    following fields in it:

    Inputs:
        S.windspeeds - vector of windspeeds
        S.all_MBC - cell array of structures with linearized OpenFAST models
                  after the MBC transfor. SS models correspond to windspeeds
        S.ws - desired windspeed for transfer function
    %}

    windspeeds  = S.windspeeds;
    all_MBC     = S.all_MBC; 
    ws          = S.ws; 
    A_set = zeros(length(windspeeds), size(all_MBC{1}.AvgA,1),  size(all_MBC{1}.AvgA,2));
    B_set = zeros(length(windspeeds), size(all_MBC{1}.AvgB,1),  size(all_MBC{1}.AvgB,2));
    C_set = zeros(length(windspeeds), size(all_MBC{1}.AvgC,1),  size(all_MBC{1}.AvgC,2));
    D_set = zeros(length(windspeeds), size(all_MBC{1}.AvgD,1),  size(all_MBC{1}.AvgD,2));
    pitch_op_set = zeros(length(windspeeds));

    for i = 1:length(windspeeds)
        A_set(i,:,:) = all_MBC{i}.AvgA;
        B_set(i,:,:) = all_MBC{i}.AvgB;
        C_set(i,:,:) = all_MBC{i}.AvgC;
        D_set(i,:,:) = all_MBC{i}.AvgD;
        pitch_op_set(i) = all_MBC{i}.ops.pitch;
    end

    A_op = interp1(windspeeds,A_set,ws);
    A_op = reshape(A_op, size(A_op,2), size(A_op,3));
    B_op = interp1(windspeeds,B_set,ws);
    B_op = reshape(B_op, size(B_op,2), size(B_op,3));
    C_op = interp1(windspeeds,C_set,ws);
    C_op = reshape(C_op, size(C_op,2), size(C_op,3));
    D_op = interp1(windspeeds,D_set,ws);
    D_op = reshape(D_op, size(D_op,2), size(D_op,3));

    ss_op = ss(A_op, B_op, C_op, D_op, ...
        'StateName', all_MBC{1}.DescStates, 'InputName', all_MBC{1}.DescInps, 'OutputName', all_MBC{1}.DescOutputs);
    
    input_string = 'collective blade-pitch';
    % input_string = 'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s'
    uind = find(contains(all_MBC{1}.DescInps, input_string));

    % output_string = 'ED PtfmPitch, (deg)';
    output_string = 'GenSpeed';
    yind = find(contains(all_MBC{1}.DescOutputs, output_string));

    Gs = tf(ss_op(yind,uind));

    fprintf('windspeed = %2.1f\n', ws)
    
end

%% Get CL transfer function
function [Hs, Gs, Cs] = get_tf_cl(S)
    %{ 
    Get closed loop transfer functions

    Inputs:
        windspeeds  - vector of windspeeds for plant 
        all_MBC     - cell array of structures with linearized OpenFAST models
                      after the MBC transfor. SS models correspond to windspeeds
        ws_tune     - vector of windspeeds for controller model 
        A_tune      - vector of "A" values from ROSCO's tuning processes
        B_tune      - vector of "B" values from ROSCO's tuning processes
        omega       - desired natural frequency
        zeta        - desired damping ratio
        ws          - desired windspeed for transfer function
    %}

    fprintf('windspeed = %2.1f\n', S.ws)
    fprintf('omega = %1.2f\n', S.omega)
    fprintf('zeta = %1.2f\n', S.zeta)

    windspeeds  = S.windspeeds;
    all_MBC     = S.all_MBC;
    ws_tune     = S.ws_tune; 
    A_tune      = S.A_tune;
    B_tune      = S.B_tune; 
    omega       = S.omega; 
    zeta        = S.zeta;
    ws          = S.ws;
    
    % --- Get set of state matrices ---
    A_set = zeros(length(windspeeds), size(all_MBC{1}.AvgA,1),  size(all_MBC{1}.AvgA,2));
    B_set = zeros(length(windspeeds), size(all_MBC{1}.AvgB,1),  size(all_MBC{1}.AvgB,2));
    C_set = zeros(length(windspeeds), size(all_MBC{1}.AvgC,1),  size(all_MBC{1}.AvgC,2));
    D_set = zeros(length(windspeeds), size(all_MBC{1}.AvgD,1),  size(all_MBC{1}.AvgD,2));
    pitch_op_set = zeros(1, length(windspeeds));

    for i = 1:length(windspeeds)
        A_set(i,:,:) = all_MBC{i}.AvgA;
        B_set(i,:,:) = all_MBC{i}.AvgB;
        C_set(i,:,:) = all_MBC{i}.AvgC;
        D_set(i,:,:) = all_MBC{i}.AvgD;
        pitch_op_set(i) = all_MBC{i}.ops.pitch;
    end

    % --- Find interpolated state space ---
    A_op = interp1(windspeeds,A_set,ws);
    A_op = reshape(A_op, size(A_op,2), size(A_op,3));
    B_op = interp1(windspeeds,B_set,ws);
    B_op = reshape(B_op, size(B_op,2), size(B_op,3));
    C_op = interp1(windspeeds,C_set,ws);
    C_op = reshape(C_op, size(C_op,2), size(C_op,3));
    D_op = interp1(windspeeds,D_set,ws);
    D_op = reshape(D_op, size(D_op,2), size(D_op,3));
    
    % Compute tower-top to platform pitch feedback
    k_float = -7.67622000000;
    select_ptfm_pitch = ismember(all_MBC{1}.DescStates, 'ED First time derivative of Platform pitch tilt rotation DOF (internal DOF index = DOF_P), rad/s')';
    k_mat = [zeros(1, size(A_op, 1));
             k_float * select_ptfm_pitch;
             zeros(1, size(A_op, 1))];
    A_op_old = A_op;
    A_op = A_op + B_op * k_mat;
    
    ss_op = ss(A_op, B_op, C_op, D_op, ...
        'StateName', all_MBC{1}.DescStates, 'InputName', all_MBC{1}.DescInps, 'OutputName', all_MBC{1}.DescOutputs);
    
    input_string = 'collective blade-pitch';
    % input_string = 'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s'
    uind = find(contains(all_MBC{1}.DescInps, input_string));

    % output_string = 'ED PtfmPitch, (deg)';
    output_string = 'GenSpeed';
    yind = find(contains(all_MBC{1}.DescOutputs, output_string));
    
    % --- System transfer function ---
    Gs = tf(ss_op(yind,uind));

    
%     % --- Get Controller "state space" ---
%     A_rosco = interp1(ws_tune, A_tune, ws);
%     B_rosco = interp1(ws_tune, B_tune, ws);
% 
%     % --- Calculate Gains --- 
%     kp_pert = 1/B_rosco * (2*zeta*omega + A_rosco);
%     ki_pert = omega^2 / B_rosco;
    
    % Get tuning values from c4 gain schedule
    blpitch_op = interp1(windspeeds,pitch_op_set,ws);
    c4_k = interp1(S.c4_oppoint', [S.c4_kp' S.c4_ki'], blpitch_op);
    kp_pert = c4_k(1);
    ki_pert = c4_k(2);

    fprintf('kp = %1.6f\n', kp_pert)
    fprintf('ki = %1.6f\n', ki_pert)
    
    % --- Controller ---
    Cs = tf([kp_pert, ki_pert],[1, 0]);
    
    % --- Get Closed Loop ---  
%     Hs = minreal(feedback(Gs*Cs, 1), .01);
    Hs = feedback(Gs*Cs, 1);
end

%% Get OL transfer function
function [Hs, Gs, Cs] = get_tf_ol(S)
    %{ 
    Get closed loop transfer functions

    Inputs:
        windspeeds  - vector of windspeeds for plant 
        all_MBC     - cell array of structures with linearized OpenFAST models
                      after the MBC transfor. SS models correspond to windspeeds
        ws_tune     - vector of windspeeds for controller model 
        A_tune      - vector of "A" values from ROSCO's tuning processes
        B_tune      - vector of "B" values from ROSCO's tuning processes
        omega       - desired natural frequency
        zeta        - desired damping ratio
        ws          - desired windspeed for transfer function
    %}

    fprintf('windspeed = %2.1f\n', S.ws)
    fprintf('omega = %1.2f\n', S.omega)
    fprintf('zeta = %1.2f\n', S.zeta)
    
    windspeeds  = S.windspeeds;
    all_MBC     = S.all_MBC;
    ws_tune     = S.ws_tune; 
    A_tune      = S.A_tune;
    B_tune      = S.B_tune; 
    omega       = S.omega; 
    zeta        = S.zeta;
    ws          = S.ws;
    
    % --- Get set of state matrices ---
    A_set = zeros(length(windspeeds), size(all_MBC{1}.AvgA,1),  size(all_MBC{1}.AvgA,2));
    B_set = zeros(length(windspeeds), size(all_MBC{1}.AvgB,1),  size(all_MBC{1}.AvgB,2));
    C_set = zeros(length(windspeeds), size(all_MBC{1}.AvgC,1),  size(all_MBC{1}.AvgC,2));
    D_set = zeros(length(windspeeds), size(all_MBC{1}.AvgD,1),  size(all_MBC{1}.AvgD,2));
    pitch_op_set = zeros(length(windspeeds));

    for i = 1:length(windspeeds)
        A_set(i,:,:) = all_MBC{i}.AvgA;
        B_set(i,:,:) = all_MBC{i}.AvgB;
        C_set(i,:,:) = all_MBC{i}.AvgC;
        D_set(i,:,:) = all_MBC{i}.AvgD;
        pitch_op_set(i) = all_MBC{i}.ops.pitch;
    end

    % --- Find interpolated state space ---
    A_op = interp1(windspeeds,A_set,ws);
    A_op = reshape(A_op, size(A_op,2), size(A_op,3));
    B_op = interp1(windspeeds,B_set,ws);
    B_op = reshape(B_op, size(B_op,2), size(B_op,3));
    C_op = interp1(windspeeds,C_set,ws);
    C_op = reshape(C_op, size(C_op,2), size(C_op,3));
    D_op = interp1(windspeeds,D_set,ws);
    D_op = reshape(D_op, size(D_op,2), size(D_op,3));
    ss_op = ss(A_op, B_op, C_op, D_op, ...
        'StateName', all_MBC{1}.DescStates, 'InputName', all_MBC{1}.DescInps, 'OutputName', all_MBC{1}.DescOutputs);
    
    input_string = 'collective blade-pitch';
    % input_string = 'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s'
    uind = find(contains(all_MBC{1}.DescInps, input_string));

    % output_string = 'ED PtfmPitch, (deg)';
    output_string = 'GenSpeed';
    yind = find(contains(all_MBC{1}.DescOutputs, output_string));
    
    % --- System transfer function ---
    Gs = tf(ss_op(yind,uind));

    
    % --- Get Controller "state space" ---
    A_rosco = interp1(ws_tune, A_tune, ws);
    B_rosco = interp1(ws_tune, B_tune, ws);

    % --- Calculate Gains --- 
    kp_pert = 1/B_rosco * (2*zeta*omega + A_rosco);
    ki_pert = omega^2 / B_rosco;

    kp_pert
    ki_pert
    
    % --- Controller ---
    Cs = tf([kp_pert, ki_pert],[1, 0]);
    
    % --- Get Open Loop ---  
    Hs = minreal(Cs*Gs, .01);
end

%% Slider functions
% Callback for all sliders defined above
function S = SliderCB(slider, eventdata, Param)
    S = guidata(slider);  % Get S structure from the figure
    S.(Param) = get(slider, 'Value');  % Any of the 'a', 'b', etc. defined
    update(S);  % Update the plot values
    guidata(slider, S);  % Store modified S in figure
end

function S = update(S)
    persistent prev_options
    [S.Hs, S.Gs, S.Cs] = S.fun(S);  
    S.h = S.pltfun(S.Hs);  % Replace old plot with new plotting values
%     S.h = S.pltfun(S.Hs, 1, 'x');  % Replace old plot with new plotting values
    hold on;
    S.pltfun(S.Hs, 1, 'x');
    hold off;
    options = getoptions(S.h);
    options.XLim = S.pltopts.XLim;
    options.YLim = S.pltopts.YLim;
%     options.XLim = prev_options.XLim;
%     options.YLim = prev_options.YLim;
    setoptions(S.h,options);
    
    % Formatting
    pfinfo = functions(S.pltfun);
    if strcmp(pfinfo.function, 'nyquistplot')
        setoptions(S.h, 'ShowFullContour', 'off');
        grid('on');
        legend(S.ax,strcat('v = ',num2str(S.ws),', \omega = ',num2str(S.omega),', \zeta = ',num2str(S.zeta)), 'fontsize',12);
    elseif strcmp(char(S.fun),'get_tf')
        legend(S.ax, strcat('v = ',num2str(S.ws)), 'fontsize',12);
    elseif strcmp(char(S.fun),'get_tf_cl')
        legend(S.ax,strcat('v = ',num2str(S.ws),', \omega = ',num2str(S.omega),', \zeta = ',num2str(S.zeta)), 'fontsize',12);
    end
end

