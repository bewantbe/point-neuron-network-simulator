% Bug: git checkout 722b6e98
% miss detect spike
addpath([getenv('HOME') '/code/point-neuron-network-simulator/mfile/']);

% bad case parameters:
%{
load('ISI_LIF-GH_ps=2mV_prps=0.16-4.5mVkHz_t=1.00e+07.mat');

s_freq = zeros(size(s_data));
s_id = 1 : numel(s_data);
for k = s_id
  s_freq(k) = 1000 ./ s_data{k}.ISI;
end

id_bad = find( diff(s_freq) > 10 );

id_bad

figure(1);
plot(s_id, s_freq, '-o', id_bad, s_freq(id_bad), 'x');

id = id_bad(3);
job = s_jobs{id};

pm = [];
%pm.prog_path = '/home/xyy/code/point-neuron-network-simulator/bin/gen_neu_dbg';
%pm.neuron_model = in_const_data.pm.neuron_model;
%pm.simu_method = in_const_data.pm.simu_method;
%pm.net     = in_const_data.pm.net;
%pm.ps_mV   = job.ps_mV;
%pm.pr      = job.pr;
%%pm.t       = 7221.5e3 * 0.5;
%pm.t       = job.t;

fprintf('fr get: %.3g Hz (spikes: %d)\n', s_freq(id), s_freq(id)/1000*job.t);

%}

pm.neuron_model = 'LIF-GH';
pm.simu_method = 'simple';
pm.net     = 'net_1_0';
pm.ps      = 3.7043000956650002e-02;  % pm.ps_mV = 2 mV
pm.pr      = 1.7374874371859297e+00;
pm.t       = 7225e3 * 0.5;
pm.seed    = [3645637114 1772554011];
pm.extra_cmd = '-v --verbose-echo';

[X, ISI, ras, pm] = gen_neu(pm, 'rm');

fprintf('fr = %.3g Hz\n', 1000 ./ ISI);

figure(2);
bg = 7221.45e3;
ids = 0:60;
plot(ids, X(1, bg+ids), '-o', 53, X(1, bg+53), 'x')

