% Using neu_psp_test

pm = [];
pm.prog_path = '../bin/gen_neu';

%s_neuron_model = {'LIF-G', 'LIF-GH'};

% should return similar ISI
s_neuron_model = {'HH-G', 'HH-GH', 'HH-GH-cont-syn'};

for id_nm = 1:length(s_neuron_model)
    pm.neuron_model = s_neuron_model{id_nm};

    PSP = neu_psp_test(pm);

    pm.simu_method  = 'auto';
    pm.net  = ones(15);
    pm.nI   = 0;
    pm.scee = 1.0 * PSP.mV_scee;  % * (PSP.t_scee/2.5);
    pm.pr   = 1.0;
    pm.ps   = 1.0 * PSP.mV_ps;
    pm.t    = 1e4;
    pm.dt   = 1.0/32;
    pm.stv  = 0.5;

    [X, ISI] = gen_neu(pm, 'new,rm');
	
    fprintf('spike freq = %.1f\n', mean(1000 ./ ISI));
end
