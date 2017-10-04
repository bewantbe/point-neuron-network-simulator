% Convergence across different simulator
addpath('../mfile');
maxabs = @(x) max(abs(x(:)));

% single neuron case

%pm0 = [];
%pm0.t    = 1e3;
%pm0.dt   = 1/32;
%pm0.stv  = pm.dt;
%pm0.pr   = 1.0;
%pm0.ps   = 0.04;
%pm0.seed = 123; % randi(2^31-1, 1, 1);

% network case

pm0 = [];
pm0.net  = ones(15);
pm0.nI   = 0;
pm0.scee = 0.05;
pm0.scie = 0.00;
pm0.scei = 0.00;
pm0.scii = 0.00;
pm0.pr   = 1.0;
pm0.ps   = 0.02;
pm0.t    = 1e3;
pm0.dt   = 1.0/1024;
pm0.stv  = pm.dt;
pm0.seed = 4563;

% legancy simulator
pm_l = pm0;
pm_l.prog_path = '../external_program/raster_tuning_HH3_gcc';
pm_l.save_poisson_events = 'tmp_poisson_events.txt';
[X_l, isi_l, ras_l, pm_expand_l] = gen_neu(pm_l, 'new,rm');
X_l(:, 1) = [];  % somehow the first data point is incorrect.

% new simulator
pm = pm0;
pm.neuron_model = 'HH-GH-cont-syn';
pm.simu_method  = 'auto';
pm.spike_threshold = 1.0;
  ie = load(pm_l.save_poisson_events);
  ie(:, 1) = ie(:, 1) + 1;
  QuietDelete(pm_l.save_poisson_events);
pm.input_event = ie;
[X, isi, ras, pm_expand] = gen_neu(pm, 'new,rm');
X(:, end) = [];  % to match X_l

pick_ras_local = @(ras, rg) ras(rg(1)*pm0.stv<ras(:,2) & ras(:,2)<rg(end)*pm0.stv, :);

rg = 1:100000;
figure(1);
ras_local = pick_ras_local(ras, rg);
ras_l_local = pick_ras_local(ras_l, rg);
plot(rg, X(:, rg),...
     rg, X_l(:, rg),...
     rg, 10*SpikeTrains(ras_local, size(X,1), length(rg), pm0.stv),...
     rg, 10*SpikeTrains(ras_l_local, size(X_l,1), length(rg), pm0.stv)...
     );

figure(2);
plot(X(rg) - X_l(rg))

disp('volt diff(rel)');
maxabs(X - X_l) / maxabs(X_l)

disp('ras diff');
maxabs(ras - ras_l)

