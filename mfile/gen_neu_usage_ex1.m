% Usage example.
% Single neuron simulation: classical HH model.

pm = [];
pm.neuron_model = 'HH-PT-GH';
pm.simu_method  = 'simple';
pm.net   = 'net_1_0';
pm.pr    = 1.0;
pm.ps_mV = 1.0;
pm.t     = 1e5;
pm.dt    = 1.0/32;
pm.stv   = pm.dt;
pm.seed  = 235478;
pm.extra_cmd = '-v --verbose-echo';

[X, ISI, ras, pm_expand, extra_data] = gen_neu(pm, 'rm,extra_data');

fprintf('Mean firing rate: %.4g Hz.\n', 1000 ./ ISI);

figure(3);
title('ISI distribution');
hist(diff(ras(ras(:,1)==1, 2)), 200);

% Show only part of data
rg = 1:floor(100/pm.stv);
s_t = rg*pm.stv;

figure(20);  % For single neuron
title('V');
ras_local = ras_pick(ras, [], s_t(end));
plot(s_t, X(:,rg)', ras_local(:,  2), 95*ones(length(ras_local),1), 'r+', 'markersize', 20);

legend('V');

figure(21);  % For single neuron
title('G');
plot(s_t, extra_data.G(:,rg)');
legend('G E', 'G I');

figure(22);  % For single neuron
title('gatings');
plot(s_t, extra_data.gatings(:,rg)');
legend('h', 'm', 'n');
% See single_neuron_dynamics.h: struct Ty_HH_GH_CUR_core for the order of gating variables (id_h, id_m, id_n)

