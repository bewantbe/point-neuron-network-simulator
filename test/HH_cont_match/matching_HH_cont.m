% s = 0.7145, 0.7278, 0.7261, 0.7369         pr=0.2*2:0.7545
% d = 0.34, 0.3371, 0.3343

function v = toMin1(d, s)

clear('pm');
pm.prog_path = '../bin/gen_neu';

pm.net  = ones(15);
pm.nI   = 0;
% 0.746/10;%delay=0.34 % 0.7145/10;% HH-GH simple  % 1.3/10; % cont-syn
%pm.scee = 0.7145/10;
pm.scie = 0.00;
pm.scei = 0.00;
pm.scii = 0.00;
pm.pr   = 0.2;
pm.ps   = 0.05;
pm.t    = 1e4;
%pm.dt   = 1.0/32;
pm.dt   = 1.0/32;
pm.stv  = 0.5;
pm.seed = 23547812;
pm.extra_cmd = '';

pm.neuron_model = 'HH-GH-cont-syn';  pm.simu_method  = 'cont-syn';  pm.scee = 1.3/10;
[X_ref, ISI_ref, ras_ref] = gen_neu(pm, 'new,rm');
X_ref = X_ref*10;

pm.neuron_model = 'HH-GH'; pm.simu_method  = 'big-delay'; pm.synaptic_delay = d;  pm.scee = s/10;
X2 = gen_neu(pm, 'new,rm');
v = norm(mean(X2) - mean(X_ref), 1);

