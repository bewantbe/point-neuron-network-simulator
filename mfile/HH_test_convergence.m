%

pm = [];
pm.prog_path = '../bin/gen_neu';

%pm.neuron_model = 'LIF-GH'; pm.simu_method  = 'simple';  pm.scee = 0.7145/20;

%pm.neuron_model = 'HH-GH'; pm.simu_method  = 'simple';  pm.scee = 0.035;

%pm.neuron_model = 'HH-GH'; pm.simu_method  = 'SSC';  pm.scee = 0.035;

pm.neuron_model = 'HH-GH'; pm.synaptic_delay = 0.2493;  pm.scee = 0.04;

%pm.neuron_model = 'HH-GH-cont-syn';  pm.simu_method  = [];  pm.scee = 0.06;

%pm.neuron_model = 'HH-PT-GH';  pm.scee = 0.04;  % ? low accuracy
%pm.simu_method  = 'SSC-Sparse';
%pm.simu_method  = 'SSC';

pm.net  = ones(15);
pm.nI   = 0;
%pm.scee = 0.04;
pm.pr   = 1.0;
pm.ps   = 0.05;
pm.t    = 1e2;
pm.dt   = 1.0/32;
pm.stv  = pm.dt;
pm.seed = 24;

pm.dt = 2 ^ -13;
[X0, ISI0, ras0] = gen_neu(pm, 'rm');

figure(2);
cla;
hd = ras_plot(ras0);
set(hd, 'linewidth', 2);
xlabel('time (ms)');
ylabel('neuron id');

figure(3);
s_t = (1:length(X0))*pm.stv;
plot(s_t, X0(1, :));

v0 = X0(:, end);
t_last = lastRASEvent(ras0, size(pm.net,1));
v0(pm.t - t_last < 4 | v0 > 10 | v0 < 0) = nan;

maxabs = @(X) max(abs(X(:)));

s_dt = 2 .^ [-6:-1:-12];
s_err = zeros(size(s_dt));
for k = 1:length(s_dt)
  pm.dt = s_dt(k);
  [X, ISI, ras] = gen_neu(pm, 'rm');

  v = X(:, end);
  t_last = lastRASEvent(ras, size(pm.net,1));
  v(pm.t - t_last < 4 | v > 10 | v < 0) = nan;
%pm.seed = 24; <4 browken, 
  %s_err(k) = maxabs(X(:, end) - X0(:, end));
  s_err(k) = maxabs(v - v0);
end

figure(1);
%loglog(s_dt, s_err, '-o');
plot(log10(s_dt), log10(s_err), '-o');

