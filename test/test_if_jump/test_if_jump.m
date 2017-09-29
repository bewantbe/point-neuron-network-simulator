%

pm = [];
pm.prog_path = '../bin/gen_neu';
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
[V, ISI, ras, pm] = gen_neu(pm, '');

%figure(1);
%plot(V');
%figure(1);
%plot(V');
