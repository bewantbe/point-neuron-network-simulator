% Demo of using configure file as parameter input.
% Extra parameters specified in "pm" will override settings in the config file.

pm = [];
pm.parameter_path = 'example_config.ini';
[X, ISI, ras] = gen_neu(pm, 'extra_data');

stv = 0.5;  % need to specify, it is not in "pm".

figure(20);
title('V');
rg = 1 : floor(100 / stv);
s_t = rg * stv;
plot(s_t, X(:,rg)');
legend('V');

