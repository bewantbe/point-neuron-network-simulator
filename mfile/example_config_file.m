% Demo of using configure file as parameter input.
% When using "parameter_path," at least pm.p needs to be provided.

pm = [];
pm.parameter_path = 'example_config.ini';
[X, ISI, ras] = gen_neu(pm, 'extra_data');

stv = 0.5;

figure(20);  % For single neuron
title('V');
rg = 1 : floor(100 / stv);
s_t = rg * stv;
plot(s_t, X(:,rg)');
legend('V');

