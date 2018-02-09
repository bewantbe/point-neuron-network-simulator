%

randMT19937('state', 123123);

pm = [];
pm.prog_path = '../../bin/gen_neu';
pm.neuron_model = 'IF-jump';
pm.simu_method = 'IF-jump';
%pm.net  = randMT19937(1000) < 0.1;
%pm.nI   = 500;
%pm.scee = 0.03;
%pm.scie = 0.02;
%pm.scei = 0.05;
%pm.scii = 0.04;
pm.net  = sparse(randMT19937(10000) < 0.1);
pm.nI   = 5000;
pm.scee = 0.005;
pm.scie = 0.008;
pm.scei = 0.009;
pm.scii = 0.010;
pm.ps   = 0.125;
pm.pr   = 2.0;
pm.t    = 1e3;
pm.dt   = 10.0;
pm.stv  = pm.dt;
pm.seed = 24;
pm.extra_cmd = '-v --verbose-echo';

tic
[V1, ISI1, ras1, pm1] = gen_neu(pm, 'new', 'data/ref_');
toc
fprintf('\n\n');
%fflush(stdout);

pm.neuron_model = 'IF-jump';
pm.simu_method = 'IF-jump-delay';
pm.synaptic_net_delay = sparse(zeros(size(pm.net)));
tic
[V2, ISI2, ras2, pm2] = gen_neu(pm, 'new');
toc
fprintf('\n\n');
%fflush(stdout);

%[ras1(18:25, :), ras2(18:25, :)]

%V1'(18:25, :) - V2'(18:25, :)

size(ras1)
size(ras2)

assert(any(size(ras1) == size(ras2)));
assert(any(ras1(:) == ras2(:)));

%figure(1);
%plot(V');

%figure(2);
%ras_plot(ras, true);

%{

find(ras1(1:1e5,1)-ras2(1:1e5,1))(1)

[ras1(72216:72229, :), ras2(72216:72229, :)]
   8770.000     12.349   8770.000     12.349
   8778.000     12.349   8778.000     12.349
   8785.000     12.349   8785.000     12.349
   8797.000     12.349   8797.000     12.349
   8819.000     12.349   8819.000     12.349
   8841.000     12.349   8820.000     12.349  <-
   8845.000     12.349   8841.000     12.349
   8854.000     12.349   8845.000     12.349
   8889.000     12.349   8854.000     12.349
   8890.000     12.349   8889.000     12.349
   8909.000     12.349   8890.000     12.349
   8923.000     12.349   8909.000     12.349
   8960.000     12.349   8923.000     12.349
   8966.000     12.349   8960.000     12.349

find(V1'(1, :)-V2'(1, :))
%}


