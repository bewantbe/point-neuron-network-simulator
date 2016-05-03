%

clear('pm');
pm.prog_path = '../bin/gen_neu';
%pm.prog_path = '/home/xyy/code/point-neuron-network-simulator-testing/bin/gen_neu';
pm.neuron_model = 'HH-GH';
pm.net  = 'net_1_0';
pm.nI   = 0;
pm.scee = 0.05;
pm.scie = 0.00;
pm.scei = 0.00;
pm.scii = 0.00;
pm.pr   = 0.05;
pm.ps   = 0.1;
pm.t    = 1e6;
pm.dt   = 1.0/32;
pm.stv  = 0.5;
pm.seed = 235478;
pm.extra_cmd = '';

tic;
[X, ISI, ras, ~, extra_data] = gen_neu(pm, 'new, extra_data');
toc;
return
ISI
figure(3);
title('ISI distribution');
hist(diff(ras(ras(:,1)==1, 2)), 200);
xlim([0, 200]);


t = (1:pm.t/pm.stv)*pm.stv;

figure(20);  % For single neuron
title('V');
plot(t, X');
legend('V');

figure(21);  % For single neuron
title('G');
plot(t, extra_data.G');
legend('G E', 'G I');

figure(22);  % For single neuron
title('gatings');
plot(t, extra_data.gatings');
legend('h', 'm', 'n');


return

old_X = X;
old_ras = ras;

maxabs = @(X) max(abs(X(:)));

maxabs(X - old_X)
maxabs(ras - old_ras)

