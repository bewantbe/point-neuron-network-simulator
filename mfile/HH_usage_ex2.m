% Usage example.
% Single neuron simulation: With sine current input.

pm = [];
pm.neuron_model = 'HH-GH-sine';
pm.simu_method  = 'simple';
pm.net   = 'net_1_0';
pm.pr    = 1.0;
pm.ps_mV = 0;
pm.t     = 1e5;
pm.dt    = 1.0/32;
pm.stv   = 0.5;
pm.seed  = 235478;
pm.extra_cmd = '-v --verbose-echo';
pm.sine_amp = 1;
pm.sine_w   = 2*pi*0.01;

[X, ISI, ras, pm] = gen_neu(pm, 'new,rm');

s_t = (1:pm.t/pm.stv)*pm.stv;

figure(20);  % For single neuron
title('V');
plot(s_t, X' - X(1));
legend('V');

