% 100 neuron test

clear('pm');
pm.prog_path = '../bin/gen_neu';

%pm.neuron_model = 'LIF-GH'; pm.simu_method  = 'simple';  pm.scee = 0.7145/20;
pm.neuron_model = 'HH-GH'; pm.simu_method  = 'simple';  pm.scee = 0.7145/10;
%pm.neuron_model = 'HH-GH'; pm.simu_method  = 'SSC';
%pm.neuron_model = 'HH-GH'; pm.simu_method  = 'SSC-Sparse';
%pm.neuron_model = 'HH-FT-GH'; pm.simu_method  = 'SSC';
%pm.neuron_model = 'HH-GH'; pm.simu_method  = 'big-delay'; pm.synaptic_delay = 0.34;  pm.scee = 0.746/10;
%pm.neuron_model = 'HH-GH-cont-syn';  pm.simu_method  = 'cont-syn';  pm.scee = 1.3/10;

% Speed
%{
simple    : 19.225009
cont-syn  : 19.648011  % ISI = 67.2132 +- 0.2
big-delay : 19.707691
SSC       : 20.110673
%}

pm.net  = ones(15);
pm.nI   = 0;
% 0.746/10;%delay=0.34 % 0.7145/10;% HH-GH simple  % 1.3/10; % cont-syn
%pm.scee = 0.7145/10;
pm.scie = 0.00;
pm.scei = 0.00;
pm.scii = 0.00;
pm.pr   = 0.2;
pm.ps   = 0.05/2;
pm.t    = 1e5;
%pm.dt   = 1.0/32;
pm.dt   = 1.0/32;
pm.stv  = 0.5;
pm.seed = 235478;
pm.extra_cmd = '';

tic;
[X, ISI, ras, ~, extra_data] = gen_neu(pm, 'new,rm');
toc;
return
if any(any(X==NaN))
  error('Wrong result!');
end

mean(ISI)
std(ISI)/sqrt(15)

figure(1)
title(sprintf('%s :: %s', pm.neuron_model, pm.simu_method));
hd = ras_plot(ras, 0, pm.t, 1:size(X,1));
set(hd, 'color', 'green');


pm.neuron_model = 'HH-GH-cont-syn';  pm.simu_method  = 'cont-syn';  pm.scee = 1.3/10;
[X_ref, ISI_ref, ras_ref] = gen_neu(pm, 'new,rm');
X_ref = X_ref*10;
%figure(4);
%t = (1:size(X_ref,2))*pm.stv;
%plot(t, mean(X_ref), t, mean(X));
norm(mean(X) - mean(X_ref), 1)

pm.neuron_model = 'HH-GH'; pm.simu_method  = 'big-delay'; pm.synaptic_delay = 0.34;  pm.scee = 0.746/10;
X2 = gen_neu(pm, 'new,rm');
norm(mean(X2) - mean(X_ref), 1)

figure(5);
t = (1:size(X_ref,2))*pm.stv;
plot(t, mean(X) - mean(X_ref), t, mean(X2) - mean(X_ref));


return

maxabs = @(X) max(abs(X(:)));

old_ras = ras;
old_X = X;

maxabs(X - old_X)
maxabs(ras - old_ras)

