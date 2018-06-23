% Convergence across for finer dt

addpath('../mfile');
maxabs = @(x) max(abs(x(:)));

% network case

pm0 = [];
pm0.neuron_model = 'LIF-GH';
pm0.simu_method  = 'SSC';

pm0.net  = ones(5);
pm0.nI   = 2;
pm0.scee_mV = 0.5;
pm0.scie_mV = 0.5;
pm0.scei_mV = 0.5;
pm0.scii_mV = 0.5;
pm0.pr      = 1.0;
pm0.ps_mV   = 1.0;
% The pm0.t should be small, otherwise the error can be too large to see
% the convergence behaviour.
pm0.t    = 1e3;
pm0.dt   = 1.0/1024;
pm0.stv  = 1.0/32;
pm0.seed = 4563;

s_dt = 1 ./ [32 64 128 256 512 1024];
s_err_V    = zeros(size(s_dt));   % error of overall volt
s_err_Vend = zeros(size(s_dt));   % error of volt at T-end
s_err_ras  = zeros(size(s_dt));   % error of last spike event timing

pm = pm0;
p = length(pm.net);

% Reference answer
pm.dt = min(s_dt)/4;
[X_r, isi_r, ras_r, pm_expand_r] = gen_neu(pm, 'new,rm');
s_se_r = arrayfun(@(id) ras_r(ras_r(:,1)==id, 2), 1:p, 'UniformOutput', false);

v0 = X_r(:, end);
t_last = lastRASEvent(ras_r, size(pm.net,1));
v0(pm.t - t_last < 4 | v0 > 10 | v0 < 0) = nan;

for id_dt = 1:numel(s_dt)
  pm.dt = s_dt(id_dt);
  [X, isi, ras, pm_expand] = gen_neu(pm, 'new,rm');

  v = X(:, end);
  t_last = lastRASEvent(ras, size(pm.net,1));
  v(pm.t - t_last < 4 | v > 10 | v < 0) = nan;

  s_err_Vend(id_dt) = maxabs(v - v0);
  s_err_V(id_dt) = maxabs(X - X_r) / maxabs(X_r);
  #lmin = min([length(ras), length(ras_r)]);
  #s_err_ras(id_dt) = maxabs(ras(1:lmin, :) - ras_r(1:lmin, :));
  
  # error in ras
  s_se = arrayfun(@(id) ras(ras(:,1)==id, 2), 1:p, 'UniformOutput', false);
  ras_err = 0;
  for id=1:p
    se   = s_se{id};
    se_r = s_se_r{id};
    lmin = min([length(se), length(se_r)]);
    ras_err = ras_err + maxabs(se(1:lmin) - se_r(1:lmin));
  end
  s_err_ras(id_dt) = ras_err;

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

