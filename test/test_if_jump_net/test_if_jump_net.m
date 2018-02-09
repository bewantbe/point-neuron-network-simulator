%

pm = [];
pm.prog_path = '../../bin/gen_neu';
pm.neuron_model = 'IF-jump';
pm.simu_method = 'IF-jump';
pm.net  = ones(5);
pm.nI   = 2;
pm.scee = 0.13;
pm.scie = 0.12;
pm.scei = 0.15;
pm.scii = 0.14;
pm.ps   = 0.125;
pm.pr   = 2.0;
pm.t    = 1e3;
pm.dt   = 1.0;
pm.stv  = pm.dt;
pm.seed = 24;
pm.extra_cmd = '-v --verbose-echo';

tic
[V1, ISI1, ras1, pm1] = gen_neu(pm, 'new', 'data/ref_');
toc
fprintf('\n\n');
fflush(stdout);

%figure(1);
%plot(V0');

%figure(2);
%ras_plot(ras0, true);

pm.neuron_model = 'IF-jump';
pm.simu_method = 'IF-jump-delay';
pm.synaptic_net_delay = zeros(size(pm.net));
tic
[V2, ISI2, ras2, pm2] = gen_neu(pm, 'new');
toc
fprintf('\n\n');
fflush(stdout);

%[ras1(18:25, :), ras2(18:25, :)]

%V1'(18:25, :) - V2'(18:25, :)

assert(any(size(ras1) == size(ras2)));

assert(any(ras1 == ras2));

%figure(1);
%plot(V');

%figure(2);
%ras_plot(ras, true);
size(ras1)
size(ras2)
