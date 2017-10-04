% test gen_neu 'legancy' mode

addpath('~/code/point-neuron-network-simulator/mfile');

clear('pm');
pm.prog_path = '/home/xyy/matcode/prj_GC_clean/HH/raster_tuning_HH3_gcc';
pm.net  = 'net_2_2';
pm.scee = 0.05;
pm.ps   = 0.04;
pm.pr   = 1.6;
pm.t    = 1e4;
pm.dt   = 1/32;
pm.stv  = 0.5;
pm.seed = 1233;

% use legancy program
[X_l, ISI_l, ras_l, pm_expand_l] = gen_neu(pm, 'new,rm');
%pm_expand.cmd_str

% use gen_neu
pm.prog_path = [];
pm.neuron_model = 'HH-GH-cont-syn';
pm.extra_cmd = '-v --verbose-echo';
[X, ISI, ras, pm_expand] = gen_neu(pm, 'new,rm');
%pm_expand.cmd_str

s_t = (1:floor(pm.t/pm.stv))*pm.stv;

figure(1);
plot(s_t, X_l.');

figure(2);
plot(s_t, X.');
