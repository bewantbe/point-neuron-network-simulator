% benchmark for neuron network.
addpath('../mfile');

warning('off', 'gen_neu:model');
warning('off', 'gen_neu:pm');

c_model = {...
'', 'IF-jump', '', 0.9
'', 'LIF-G', 'SSC', 0.9
'', 'LIF-GH', 'SSC', 0.9
'', 'HH-G', 'SSC', 1.73
'', 'HH-GH', 'SSC', 1.73
'', 'HH-PT-GH', 'SSC', 1.73
'', 'HH-GH-sine', 'SSC', 1.73
'', 'HH-GH-cont-syn', '', 1.73
'', '', '', 0
'', 'LIF-G', 'simple', 0.9
'', 'LIF-GH', 'simple', 0.9
'', 'HH-G', 'simple', 1.73
'', 'HH-GH', 'simple', 1.73
'', 'HH-PT-GH', 'simple', 1.73
'', 'HH-GH-sine', 'simple', 1.73
'','','', 0
'', 'LIF-G', 'SSC-Sparse', 0.9
'', 'LIF-GH', 'SSC-Sparse', 0.9
'', 'HH-G', 'SSC-Sparse', 1.73
'', 'HH-GH', 'SSC-Sparse', 1.73
'', 'HH-PT-GH', 'SSC-Sparse', 1.73
'', 'HH-GH-sine', 'SSC-Sparse', 1.73
'','','', 0
'', 'LIF-G', 'SSC-Sparse2', 0.9
'', 'LIF-GH', 'SSC-Sparse2', 0.9
'', 'HH-G', 'SSC-Sparse2', 1.73
'', 'HH-GH', 'SSC-Sparse2', 1.73
'', 'HH-PT-GH', 'SSC-Sparse2', 1.73
'', 'HH-GH-sine', 'SSC-Sparse2', 1.73
'','','', 0
'', 'LIF-G', 'big-delay', 0.9
'', 'LIF-GH', 'big-delay', 0.9
'', 'HH-G', 'big-delay', 1.73
'', 'HH-GH', 'big-delay', 1.73
'', 'HH-PT-GH', 'big-delay', 1.73
'', 'HH-GH-sine', 'big-delay', 1.73
'','','', 0
'../external_program/raster_tuning_LIF_icc','legancy-LIF-G','', 0.9
'../external_program/raster_tuning_LIF_GH_icc','legancy-LIF-GH','', 0.9
'../external_program/raster_tuning_HH3_gcc', 'legancy-HH-GH-cont-syn', '', 1.73
};
s_model = cell2struct(c_model, {'exe_path', 'model', 'simu_method', 'prps_mV'}, 2);

randMT19937('state', 1234);

pm = [];
pm.net  = randMT19937(1000) < net_sparsity;
pm.nE   = 800;
pm.nI   = 200;
pm.t    = 1e3;
pm.dt   = 1/32;
pm.stv  = 0.5;
pm.pr      = 2.0;
pm.ps_mV   = 0.9;
pm.scee_mV = 0.004;  % Small value so that network dynamics will not change.
pm.scie_mV = 0.005;
pm.scei_mV = 0.006;
pm.scii_mV = 0.003;
pm.seed = 123;
pm_description = sprintf('Network: sparsity = %.3g',...
  sum(pm.net(:)) / (numel(pm.net) - pm.nE - pm.nI));

fprintf('Description: %s\n', pm_description);
fprintf('  n = %d+%d, t = %.3g s, dt = 1/%.3g ms, stv = %.3g ms\n',...
        pm.nE, pm.nI, pm.t/1000, 1/pm.dt, pm.stv);
fprintf('model                     sec     mean freq (Hz)');
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
    fprintf('%-26s', [pm.neuron_model ' + ' pm.simu_method]);
  else
    fprintf('%-26s', pm.neuron_model);
  end
  if length([pm.prog_path pm.neuron_model pm.simu_method])==0
    fprintf('\n');
    continue
  end
  pm.ps_mV = s_model(id_model).prps_mV / pm.pr;
  t_start = tic;
  [X, isi, ras, pm_expand] = gen_neu(pm, 'new,rm');
  fprintf('%-8.3f', toc(t_start));
  fprintf('%-8.3f', mean(1000 ./ isi));
  fprintf('\n');
end

