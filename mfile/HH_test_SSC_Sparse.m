% test SSC and SSC-Sparse
tocs = @(st) fprintf('%s: t = %6.3f s\n', st, toc());

pm = [];
pm.prog_path = '../bin/gen_neu';
%pm.neuron_model = 'HH-PT-GH';
pm.neuron_model = 'HH-GH';

randMT19937('state', [1 342 232]);
pm.net  = 1*(randMT19937(1000)<0.1);
pm.nI   = 0;
pm.scee = 0.02;
pm.pr   = 1.0;
pm.ps   = 0.03;
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


pm.simu_method  = 'SSC-Sparse';
tic
[X1, ISI1, ras1] = gen_neu(pm, 'new,rm');
toc

figure(2);
hd = ras_plot(ras0);
set(hd, 'linewidth', 2);
xlabel('time (ms)');
ylabel('neuron id');

figure(3);
s_t = (1:length(X0))*pm.stv;
plot(s_t, X0(1, :));

semilogy(sort(diff(ras0(:,2))))

figure(4)
loglog(sort(diff(ras0(:,2))))

