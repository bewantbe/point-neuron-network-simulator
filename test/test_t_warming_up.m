% Test:
% option: --t-warming-up
addpath('../mfile');

pm = [];
pm.neuron_model = 'HH-GH';
pm.simu_method  = 'auto';
pm.net   = ones(100);
pm.pr    = 1.0;
pm.ps_mV = 1.0;
pm.scee_mV = 0.3;
pm.t     = 1e3;
pm.dt    = 1.0/32;
pm.stv   = 0.5;
pm.seed  = 25478;
pm.extra_cmd = '-v --verbose-echo';

[o_X, o_ISI, o_ras, ~, o_extra_data] = gen_neu(pm, 'new,rm,extra_data,ext_T');

pm.extra_cmd = '-v --verbose-echo --t-warming-up 1e2';
[X, ISI, ras, ~, extra_data] = gen_neu(pm, 'new,rm,extra_data');

maxabs = @(x) max(abs(x(:)));

% should output 0
maxabs(X - o_X)
maxabs(extra_data.G - o_extra_data.G)
maxabs(extra_data.gatings - o_extra_data.gatings)
maxabs(ras - o_ras)

% o_ISI and ISI is expected to be different, due to defect of gen_neu.m
maxabs(o_ISI - ISI)
