% in /home/yanyang/matcode/point-neuron-network-simulator/mfile

pm = [];
pm.neuron_model = 'LIF-GH';  % one of LIF-G, LIF-GH, HH-GH
pm.nE   = 1000;
pm.nI   = 0;
pm.net  = rand(pm.nE+pm.nI)>0.8;
pm.scee = 0.001;
pm.scie = 0.00;
pm.scei = 0.00;
pm.scii = 0.00;
pm.pr   = 1.6;
pm.ps   = 0.02;
pm.t    = 200;
pm.dt   = 1.0/32;
pm.stv  = 0.5;
pm.seed = 'auto';
pm.extra_cmd = '--input-event-path ep.txt';

p = pm.nE+pm.nI;
poisson_rate = 1*p;  % kHz
elen = pm.t*poisson_rate;
a = exprnd(1/poisson_rate, elen/2, 1);
a = [a; exprnd(1/(poisson_rate/2), elen/2, 1)];
a = cumsum(a);

tic
a = [randi([0 p-1], size(a)) a];
fd = fopen('ep.txt','w');
fprintf(fd, '%d %.6f\n', a');
fclose(fd);
toc

%{
tic
fd = fopen('ep.txt','w');
for k=1:size(a,1)
  fprintf(fd, '%d %.6f\n', randi([0 p-1]), a(k));
end
fclose(fd);
toc
%}

tic
[X, ISI, ras, pm] = gen_neu(pm, 'new,rm');
toc

% figure(1);
% plot(X(1,:));

figure(1);
hist(ras(:, 2), 100);
axis([0 200 0 400]);
%xlim([0 200]);

