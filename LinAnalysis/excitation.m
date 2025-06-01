%% Parameterize Excitation Signals

%% Step
step.time = 1;
step.init_val = 0;
step.final_val = 1;
step.sample_time = 0.1;

%% Ramp
ramp.slope = 1;
ramp.start_time = 0;
ramp.init_output = 1;

%% Sine
sine.amp = 1;
sine.bias = 0;
sine.freq = 1;
sine.phase = 0;
sine.sample_time = 0.001;

%% Chirp
chirp.init_freq = 50;
chirp.target_time = 1;
chirp.target_time_freq = 100;


%% Gaussian Noise
gaussian.amp = 0.05;
gaussian.freq = 0.001;
guassian.phase = 0;