% test SSC and SSC-Sparse
tocs = @(st) fprintf('%s: t = %6.3f s\n', st, toc());

pm = [];
pm.prog_path = '../bin/gen_neu';
pm.neuron_model = 'LIF-GH';

randMT19937('state', [1 342 232]);
pm.net  = 1*(randMT19937(1000)<0.1);
pm.nI   = 0;
pm.scee = 0.0016;
pm.pr   = 1.0;
pm.ps   = 0.0125;
pm.t    = 1e3;
pm.dt   = 1.0/32;
pm.stv  = pm.dt;
pm.seed = 24;

s_simu_method  = {'simple', 'big-delay', 'SSC', 'SSC-Sparse', 'SSC-Sparse2'};
%s_simu_method  = {'simple'};

for id_sm = 1:length(s_simu_method)
    pm.simu_method = s_simu_method{id_sm};
    if strcmp(pm.simu_method, 'big-delay')
        pm.synaptic_delay = 1.01*pm.dt;
    else
        pm.synaptic_delay = [];
    end
    tic
    [X0, ISI0, ras0] = gen_neu(pm, 'new,rm');
    tocs(pm.simu_method)

    fprintf('Mean freq = %.4f\n', mean(1000 ./ ISI0));
end

return
