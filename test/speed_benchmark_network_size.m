% Network change size
addpath('../mfile');

warning('off', 'gen_neu:model');
warning('off', 'gen_neu:pm');

c_model = {...
'', 'IF-jump', '', 0.9
'', 'LIF-GH', 'SSC', 0.9
'', 'HH-GH', 'SSC', 1.73
'', '', '', 0
'', 'LIF-GH', 'simple', 0.9
'', 'HH-GH', 'simple', 1.73
'','','', 0
'', 'LIF-GH', 'SSC-Sparse', 0.9
'', 'HH-GH', 'SSC-Sparse', 1.73
'','','', 0
'', 'LIF-GH', 'SSC-Sparse2', 0.9
'', 'HH-GH', 'SSC-Sparse2', 1.73
'','','', 0
'', 'LIF-GH', 'big-delay', 0.9
'', 'HH-GH', 'big-delay', 1.73
};
s_model = cell2struct(c_model, {'exe_path', 'model', 'simu_method', 'prps_mV'}, 2);

pm = [];
pm.t    = 1e3;
pm.dt   = 1/32;
pm.stv  = 0.5;
pm.pr      = 2.0;
pm.ps_mV   = 0.9;
pm.scee_mV_1000 = 0.004;
pm.scie_mV_1000 = 0.005;
pm.scei_mV_1000 = 0.006;
pm.scii_mV_1000 = 0.003;
pm.seed = 123;
net_sparsity = 0.1;
pm_description = 'Network of changing size';

s_nn = [100 300 1000 3000 10000];

fprintf('Description: %s\n', pm_description);
fprintf('  t = %.3g s, dt = 1/%.3g ms, stv = %.3g ms\n',...
        pm.t/1000, 1/pm.dt, pm.stv);
fprintf('                        sec     mean freq (Hz)');
fprintf('\nmodel \\ nn              ');
for id_nn = 1:numel(s_nn)
  fprintf('%-16d', s_nn(id_nn));
end
fprintf('\n');
for id_model = 1:numel(s_model)
  pm.prog_path    = s_model(id_model).exe_path;
  pm.neuron_model = s_model(id_model).model;
  pm.simu_method  = s_model(id_model).simu_method;
  if strfind(pm.simu_method, 'delay')
    pm.synaptic_delay = 1.01*pm.dt;
  else
    pm.synaptic_delay = [];
  end
  if ~isempty(pm.simu_method)
    fprintf('%-24s', [pm.neuron_model ' + ' pm.simu_method]);
  else
    fprintf('%-24s', pm.neuron_model);
  end
  if length([pm.prog_path pm.neuron_model pm.simu_method])==0
    fprintf('\n');
    continue
  end
  randMT19937('state', 1234);
  pm.ps_mV = s_model(id_model).prps_mV / pm.pr;
  for id_nn = 1:numel(s_nn)
    pm.net  = sparse(randMT19937(s_nn(id_nn)) < net_sparsity);
    pm.scee_mV = pm.scee_mV_1000 * 1000/s_nn(id_nn);
    pm.scie_mV = pm.scie_mV_1000 * 1000/s_nn(id_nn);
    pm.scei_mV = pm.scei_mV_1000 * 1000/s_nn(id_nn);
    pm.scii_mV = pm.scii_mV_1000 * 1000/s_nn(id_nn);
    t_start = tic;
    [X, isi, ras, pm_expand] = gen_neu(pm, 'new,rm');
    fprintf('%-8.3f', toc(t_start));
    fprintf('%-8.3f', mean(1000 ./ isi));
  end
  
  fprintf('\n');
end

