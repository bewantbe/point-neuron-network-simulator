% single neuron test

clear('pm');
pm.prog_path = '../bin/gen_neu';
%pm.prog_path = '/home/xyy/code/point-neuron-network-simulator-testing/bin/gen_neu';
pm.simu_method  = 'simple';

pm.net  = [1];
pm.scee = 0;
pm.pr   = 2;
pm.ps   = 0.01;
pm.t    = 1e3;
pm.dt   = 1.0/16;
pm.stv  = pm.dt;
pm.seed = 235478;

t = (1:floor(pm.t/pm.dt))*pm.stv;

close all

pm.neuron_model = 'HH-G';
[V, ISI, ras, ~, extra_data] = gen_neu(pm, 'new,rm, extra_data');
figure();
plot(t, V);  title(sprintf('%s :: %s', pm.neuron_model, pm.simu_method));

pm.neuron_model = 'HH-GH';
pm.ps = pm.ps * 2;
[V, ISI, ras, ~, extra_data] = gen_neu(pm, 'new,rm, extra_data');
figure();
plot(t, V);  title(sprintf('%s :: %s', pm.neuron_model, pm.simu_method));

pm.neuron_model = 'HH-G-sine';
pm.sine_amp = 1;
pm.sine_freq   = 0.01;
pm.pr = 1.0;
pm.ps = 0.01;
[V, ISI, ras, ~, extra_data] = gen_neu(pm, 'new,rm, extra_data');
figure();
plot(t, V);  title(sprintf('%s :: %s', pm.neuron_model, pm.simu_method));

pm.neuron_model = 'HH-GH-sine';
pm.ps = pm.ps * 2;
[V, ISI, ras, ~, extra_data] = gen_neu(pm, 'new,rm, extra_data');
figure();
plot(t, V);  title(sprintf('%s :: %s', pm.neuron_model, pm.simu_method));

pm.neuron_model = 'HH-G-extI';
pm.extI = 6.5;
pm.pr = 0.5;
pm.ps = 0.01;
[V, ISI, ras, ~, extra_data] = gen_neu(pm, 'new,rm, extra_data');
figure();
plot(t, V);  title(sprintf('%s :: %s', pm.neuron_model, pm.simu_method));

pm.neuron_model = 'HH-GH-extI';
pm.ps = pm.ps * 2;
[V, ISI, ras, ~, extra_data] = gen_neu(pm, 'new,rm, extra_data');
figure();
plot(t, V);  title(sprintf('%s :: %s', pm.neuron_model, pm.simu_method));

%ras_plot(ras, 0, pm.t, 1:100);
%any(any(X==NaN))

