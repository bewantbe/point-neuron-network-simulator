% test HH model

addpath('/home/xyy/code/point-neuron-network-simulator/mfile/');

% Single neuron case

% Get reference answer
pm = [];
pm.prog_path = '/home/xyy/code/point-neuron-network-simulator/bin/gen_neu_ref';
pm.neuron_model = 'HH-PT-GH';
pm.net  = 'net_1_0';
%pm.net = ones(15);
pm.nI   = 0;
pm.scee = 0.05;
pm.pr   = 1.0;
pm.ps   = 0.04;
pm.t    = 1e3;
pm.dt   = 1.0/1024;
pm.stv  = pm.dt;
pm.seed = 4563;
tic
[X_ref, ISI_ref, ras_ref] = gen_neu(pm, 'new,rm');
toc

%
pm.prog_path = '/home/xyy/code/point-neuron-network-simulator/bin/gen_neu';
pm.neuron_model = 'HH-GH';
pm.synaptic_delay = 0.2493;

tic
[X, ISI, ras] = gen_neu(pm, 'new,rm');
toc

%maxabs = @(X) max(abs(X(:)));
%maxabs(X - X_ref) / maxabs(X_ref)
%maxabs(ras - ras_ref)

figure(5);
plot(ras(:,2) - ras_ref(:,2), 'o');
std_err = std(ras(2:end,2) - ras_ref(2:end,2))
mean_err = mean(ras(2:end,2) - ras_ref(2:end,2))

% for threshold 65mV
%   the spiking time shift (from peak) is 0.2493 +- 0.0020
% for threshold 73.35mV
%   the spiking time shift (from peak) is 0.2133 +- 0.0003 for weak input(ISI 25Hz), and 0.2128 +- 0.0022
