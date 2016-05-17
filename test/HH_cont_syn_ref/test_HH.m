% test HH model

run('/home/xyy/matcode/GC_clean/startup.m');
addpath('/home/xyy/matcode/prj_GC_clean/HH/');

addpath('/home/xyy/code/point-neuron-network-simulator/mfile/');

% Single neuron case

% Get reference answer
clear('pm');
pm.neuron_model = '/home/xyy/code/point-neuron-network-simulator/test/HH_cont_syn_ref/IFsimu_HH/bin/Release/raster_tuning_HH3_gcc';
%pm.net  = 'net_1_0';
pm.net  = ones(15);
pm.nI   = 0;
pm.scee = 0.05;
pm.scie = 0.00;
pm.scei = 0.00;
pm.scii = 0.00;
pm.pr   = 1.0;
pm.ps   = 0.02;
pm.t    = 1e3;
pm.dt   = 1.0/1024;
pm.stv  = pm.dt;
pm.seed = 4563;
pm.extra_cmd = '--save-poisson-events poisson_events.txt -q';
tic
[X_ref, ISI_ref, ras_ref] = gen_HH(pm, 'rm');
toc
X_ref(:, 1) = [];  % somehow the first data point is incorrect.

%
pm.prog_path = '../../bin/gen_neu';
pm.neuron_model = 'HH-GH-cont-syn';
pm.simu_method  = 'cont-syn';
pm.extra_cmd = '--initial-state-path init.txt --input-event-path poisson_events.txt';

tic
[X, ISI, ras] = gen_neu(pm, 'rm');
toc
X(:, end) = [];

figure(3);
id_n = 3;
rg = 1:round(pm.t/pm.stv)-1;
%rg = rg(round(end/2):end);
t = rg*pm.stv;
plot(t, X_ref(id_n,rg), '-x', t, X(id_n,rg), '-x');

figure(4);
plot(t, X_ref(id_n,rg) - X(id_n,rg), '-x');


maxabs = @(X) max(abs(X(:)));
maxabs(X - X_ref) / maxabs(X_ref)
maxabs(ras - ras_ref)

