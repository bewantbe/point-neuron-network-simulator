% Convergence across different simulator
addpath('../mfile');
maxabs = @(x) max(abs(x(:)));

c_model = {...
'../external_program/raster_tuning_HH3_gcc', 'HH-GH-cont-syn'
'../external_program/raster_tuning_HH3_gcc49', 'HH-GH-cont-syn'
'../external_program/raster_tuning_HH3_gcc49_sparse', 'HH-GH-cont-syn'
};
s_model = cell2struct(c_model, {'legancy_path', 'model'}, 2);

% single neuron case

pm0 = [];
pm0.t    = 1e3;
pm0.dt   = 1/32;
pm0.stv  = pm0.dt;
pm0.pr   = 1.0;
pm0.ps   = 0.04;
pm0.seed = 123; % randi(2^31-1, 1, 1);

% network case

%pm0 = [];
%pm0.net  = ones(15);
%pm0.nI   = 0;
%pm0.scee = 0.05;
%pm0.scie = 0.00;
%pm0.scei = 0.00;
%pm0.scii = 0.00;
%pm0.pr   = 1.0;
%pm0.ps   = 0.02;
%pm0.t    = 1e3;
%pm0.dt   = 1.0/1024;
%pm0.stv  = pm0.dt;
%pm0.seed = 4563;

s_err_V   = zeros(size(s_model));
s_err_ras = zeros(size(s_model));

for id_model = 1:numel(s_model)
  % legancy simulator
  pm_l = pm0;
  pm_l.prog_path = s_model(id_model).legancy_path;
  pm_l.save_poisson_events = 'tmp_poisson_events.txt';
  pm_l.extra_cmd = '-q';
  [X_l, isi_l, ras_l, pm_expand_l] = gen_neu(pm_l, 'new,rm');
  X_l(:, 1) = [];  % somehow the first data point is incorrect.

  % new simulator
  pm = pm0;
  pm.neuron_model = s_model(id_model).model;
  pm.simu_method  = 'auto';
  pm.spike_threshold = 1.0;
    ie = load(pm_l.save_poisson_events);
    ie(:, 1) = ie(:, 1) + 1;
    QuietDelete(pm_l.save_poisson_events);
  pm.input_event = ie;
  pm.extra_cmd = '';
  [X, isi, ras, pm_expand] = gen_neu(pm, 'new,rm');
  X(:, end) = [];  % to match X_l

  s_err_V(id_model) = maxabs(X - X_l) / maxabs(X_l);
  s_err_ras(id_model) = maxabs(ras - ras_l);
end

figure(1);
semilogy(1:numel(s_model), s_err_V, '-o');
ylabel('V err');
xlabel('dt (ms)');

figure(2);
semilogy(1:numel(s_model), s_err_ras, '-o');
ylabel('ras err');
xlabel('dt (ms)');

return

pick_ras_local = @(ras, rg) ras(rg(1)*pm0.stv<ras(:,2) & ras(:,2)<rg(end)*pm0.stv, :);

rg = 1:100000;
figure(21);
ras_local = pick_ras_local(ras, rg);
ras_l_local = pick_ras_local(ras_l, rg);
plot(rg, X(:, rg),...
     rg, X_l(:, rg),...
     rg, 10*SpikeTrains(ras_local, size(X,1), length(rg), pm0.stv),...
     rg, 10*SpikeTrains(ras_l_local, size(X_l,1), length(rg), pm0.stv)...
     );

figure(21);
plot(X(rg) - X_l(rg))
