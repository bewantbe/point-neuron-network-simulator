% Usage example.
% Single neuron simulation: classical HH model.

pm = [];
pm.neuron_model = 'HH-GH';
pm.simu_method  = 'simple';
pm.net   = 'net_1_0';
pm.pr    = 1.0;
pm.ps_mV = 1.0;
pm.t     = 1e5;
pm.dt    = 1.0/32;
pm.stv   = 0.5;
pm.seed  = 235478;
pm.extra_cmd = '-v --verbose-echo';

[X, ISI, ras, pm_expand, extra_data] = gen_neu(pm, 'rm,extra_data');

fprintf('Mean firing rate: %.4g Hz.\n', 1000 ./ ISI);

figure(3);
title('ISI distribution');
hist(diff(ras(ras(:,1)==1, 2)), 200);

s_t = (1:pm.t/pm.stv)*pm.stv;

figure(20);  % For single neuron
title('V');
plot(s_t, X');
legend('V');

figure(21);  % For single neuron
title('G');
plot(s_t, extra_data.G');
legend('G E', 'G I');

figure(22);  % For single neuron
title('gatings');
plot(s_t, extra_data.gatings');
legend('h', 'm', 'n');
% See single_neuron_dynamics.h: struct Ty_HH_GH_CUR_core for the order of gating variables (id_h, id_m, id_n)

