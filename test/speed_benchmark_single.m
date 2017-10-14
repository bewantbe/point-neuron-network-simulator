% Benchmark of single neuron simulation.
addpath('../mfile');

warning('off', 'gen_neu:model');

c_model = {...
'', 'Hawkes-GH', 'simple', 0.9
'', 'IF-jump', '', 0.9
'', 'HH-GH-cont-syn', '', 4.0
'', '', '', 0
'', 'LIF-G', 'SSC', 0.9
'', 'LIF-GH', 'SSC', 0.9
'', 'HH-G', 'SSC', 4.0
'', 'HH-GH', 'SSC', 4.0
'', 'HH-PT-GH', 'SSC', 4.0
'', 'HH-GH-sine', 'SSC', 4.0
'', '', '', 0
'', 'LIF-G', 'simple', 0.9
'', 'LIF-GH', 'simple', 0.9
'', 'HH-G', 'simple', 4.0
'', 'HH-GH', 'simple', 4.0
'', 'HH-PT-GH', 'simple', 4.0
'', 'HH-GH-sine', 'simple', 4.0
'', '', '', 0
'', 'LIF-G', 'SSC-Sparse', 0.9
'', 'LIF-GH', 'SSC-Sparse', 0.9
'', 'HH-G', 'SSC-Sparse', 4.0
'', 'HH-GH', 'SSC-Sparse', 4.0
'', 'HH-PT-GH', 'SSC-Sparse', 4.0
'', 'HH-GH-sine', 'SSC-Sparse', 4.0
'', '', '', 0
'', 'LIF-G', 'SSC-Sparse2', 0.9
'', 'LIF-GH', 'SSC-Sparse2', 0.9
'', 'HH-G', 'SSC-Sparse2', 4.0
'', 'HH-GH', 'SSC-Sparse2', 4.0
'', 'HH-PT-GH', 'SSC-Sparse2', 4.0
'', 'HH-GH-sine', 'SSC-Sparse2', 4.0
'', '', '', 0
'', 'LIF-G', 'big-delay', 0.9
'', 'LIF-GH', 'big-delay', 0.9
'', 'HH-G', 'big-delay', 4.0
'', 'HH-GH', 'big-delay', 4.0
'', 'HH-PT-GH', 'big-delay', 4.0
'', 'HH-GH-sine', 'big-delay', 4.0
'','','', 0
'../external_program/raster_tuning_LIF_icc','legancy-LIF-G','', 0.9
'../external_program/raster_tuning_LIF_GH_icc','legancy-LIF-GH','', 0.9
'../external_program/raster_tuning_HH3_gcc', 'legancy-HH-GH-cont-syn', '', 4.0
};
s_model = cell2struct(c_model, {'exe_path', 'model', 'simu_method', 'prps_mV'}, 2);

%% Single neuron test.
pm = [];
pm.net  = 1;
pm.nE   = 1;
pm.nI   = 0;
pm.t    = 1e5;
pm.dt   = 1/32;
pm.stv  = 0.5;
pm.seed = 123;
pm_description = 'Single neuron test across models and input poisson rates.';

s_pr = [1 8 64 512];

fprintf('Description: %s\n', pm_description);
fprintf('  n = %d+%d, t = %.3g s, dt = 1/%.3g ms, stv = %.3g ms\n',...
        pm.nE, pm.nI, pm.t/1000, 1/pm.dt, pm.stv);
fprintf('%-24s', 'model \ sec \ pr');
for id_pr = 1:numel(s_pr)
  fprintf('%-8.3g', s_pr(id_pr));
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
  if strfind(pm.neuron_model, 'Hawkes')
  disp('gaga');
    pm.spike_threshold = 0.03;
  else
    pm.spike_threshold = [];
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
  s_fq = zeros(size(s_pr));
  for id_pr = 1:numel(s_pr)
    pm.pr    = s_pr(id_pr);
    pm.ps_mV = s_model(id_model).prps_mV / pm.pr;  % "0.9 / pr" so that firing rate is about 31 Hz

    t_start = tic;
    [X, isi, ras, pm_expand] = gen_neu(pm, 'new,rm');
    fprintf('%-8.3f', toc(t_start));
    s_fq(id_pr) = 1000 ./ isi;
  end
  fprintf('    ');
  for fq = s_fq
    fprintf('%-8.1f', fq);
  end
  fprintf('\n');
end

