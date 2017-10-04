% Convergence across different simulator
addpath('../mfile');
maxabs = @(x) max(abs(x(:)));

% single neuron case

%pm = [];
%pm.t    = 1e3;
%pm.dt   = 1/32;
%pm.stv  = 1/32;
%pm.pr   = 1.0;
%pm.ps   = 0.04;
%pm.seed = 123; % randi(2^31-1, 1, 1);
%pm.neuron_model = 'HH-GH-cont-syn';
%pm.simu_method  = 'auto';
%pm.spike_threshold = 1.0;
%pm.input_event = ie;
%pm.extra_cmd = '';

% network case

pm0 = [];
pm0.net  = ones(15);
pm0.nI   = 0;
pm0.scee = 0.05;
pm0.scie = 0.06;
pm0.scei = 0.07;
pm0.scii = 0.04;
pm0.pr   = 1.0;
pm0.ps   = 0.02;
pm0.t    = 1e3;
pm0.dt   = 1.0/1024;
pm0.stv  = pm0.dt;
pm0.seed = 4563;

s_dt = 1 ./ [32 64 128 256 512 1024];
s_err_V   = zeros(size(s_dt));
s_err_ras = zeros(size(s_dt));

% Reference answer
pm.dt = s_dt(end)/2;
[X_r, isi_r, ras_r, pm_expand_r] = gen_neu(pm, 'new,rm');

for id_dt = 1:numel(s_dt)
  pm.dt = s_dt(id_dt);
  [X, isi, ras, pm_expand] = gen_neu(pm, 'new,rm');

  s_err_V(id_dt) = maxabs(X - X_r) / maxabs(X_r);
  s_err_ras(id_dt) = maxabs(ras - ras_r);
end

figure(1);
loglog(s_dt, s_err_V, '-o');
ylabel('V err');
xlabel('dt (ms)');

figure(2);
loglog(s_dt, s_err_ras, '-o');
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


