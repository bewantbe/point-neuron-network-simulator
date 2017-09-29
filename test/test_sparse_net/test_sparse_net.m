addpath('../mfile');

adj = rand(10) < 0.2;

pm = [];
pm.neuron_model = 'HH-GH';
pm.simu_method = 'simple';
pm.net  = sparse(adj);
pm.nI   = 0;
pm.scee_mV = 1.5;
pm.scie_mV = 0.0;
pm.scei_mV = 0.0;
pm.scii_mV = 0.0;
pm.pr      = 2.0 * ones(1,length(pm.net));
pm.ps_mV   = 1.0;
pm.t    = 1e3;
pm.dt   = 1/32;
pm.stv  = 0.5;
pm.seed = 1234;
pm.extra_cmd = '-v --verbose-echo';

[X, ISI, ras, pm1] = gen_neu(pm, 'new,rm');

% Stress test
n = 1e5;
ii = [2:n 1]';
jj = (1:n)';
val = ones(n, 1);
a = spconvert([ii jj val]);

pm = [];
pm.neuron_model = 'LIF-G';
pm.simu_method = 'simple';
pm.net  = a;
pm.nI   = 0;
pm.scee_mV = 1.0;
pm.scie_mV = 0.0;
pm.scei_mV = 0.0;
pm.scii_mV = 0.0;
pm.pr      = 1.0;
pm.ps_mV   = 1.0;
pm.t    = 1e3;
pm.dt   = 1/32;
pm.stv  = 0.5;
pm.seed = 1234;
pm.extra_cmd = '-v --verbose-echo';

[X, ISI, ras, pm1] = gen_neu(pm, 'new,rm');

