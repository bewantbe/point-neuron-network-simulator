% 100 neuron test

clear('pm');
pm.prog_path = '../bin/gen_neu';

pm.neuron_model = 'HH-GH-cont-syn';  pm.simu_method  = 'cont-syn';

pm.net  = 'net_2_2';
pm.nI   = 0;
pm.scee = 0.07;
pm.scie = 0.00;
pm.scei = 0.00;
pm.scii = 0.00;
pm.pr   = 1.0;
pm.ps   = 0.02;
pm.pr_mul = [1 0];
pm.t    = 1e4;
pm.dt   = 1.0/32/4;
pm.stv  = pm.dt;
pm.seed = 45665;
pm.extra_cmd = '';

[X_ref, ISI_ref, ras_ref] = gen_neu(pm, 'new,rm');
X_ref = X_ref*10;

%figure(1);
%t = pm.stv*(1:size(X_ref,2));
%plot(t, X_ref');

t_range = [-5 20];
id_from = 1;
id_to = 2;
[Y1, T1] = SpikeSlicing(X_ref, ras_ref, pm.stv, id_to, id_from, t_range);
%figure(11);
%plot(T1', Y1');

if any(any(X_ref==NaN))
  error('Wrong result!');
end

pm.neuron_model = 'HH-GH'; pm.simu_method  = 'big-delay'; pm.synaptic_delay = 0.34;
pm.scee = pm.scee * 0.54;
[X, ISI, ras, ~, ex] = gen_neu(pm, 'new,rm');

[Y2, T2] = SpikeSlicing(X, ras, pm.stv, id_to, id_from, t_range);
%figure(12);
%plot(T2', Y2');

figure(13);
plot(T1', Y1', 'b', T2', Y2', 'r', 'linewidth', 2);


%%figure(10);
%%plot(t, X_ref(1,:), t, X(1,:));

%figure(11);
%plot(t, X_ref(1,:) - X(1,:));

