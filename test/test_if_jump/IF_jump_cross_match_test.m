% Test consistency across versions
fprintf('~~~~~~ IF-jump Simulation Matching Test ~~~~~~\n');

addpath('../../mfile');
maxabs = @(x) max(abs(x(:)));

ref_prog_path = '~/tmp/IF_Jump/';

pm = [];
pm.neuron_model = 'IF-jump';
pm.simu_method  = 'IF-jump';
pm.net     = read_compact_sparse_net([ref_prog_path '@@test-network.txt']);
pm.nI      = length(pm.net) / 2;
pm.scee = 0.1;
pm.scie = 0.09;
pm.scei = 0.15;
pm.scii = 0.12;
pm.pr   = 3.0;
pm.ps   = 0.04;
pm.t    = 1000.0;
pm.dt   = 1.0;
pm.stv  = pm.dt;
pm.seed = 0;
pm.extra_cmd = ['--verbose-echo --input-event-path ' ref_prog_path 'poi.txt'];
[X_poi, ISI_poi, ras_poi, pm] = gen_neu(pm, 'new,rm');


% Read ref ans
ras_ref = load([ref_prog_path '@@test-13-9_8-2017-14_9_36-fire_time.txt']);
ras_ref = ras_ref(:, [2 1]);
ras_ref(:, 1) += 1;

% sort
%[~, ids] = sort();

poi = load([ref_prog_path 'poi.txt']);

%poi_local = poi(1:find(poi(:,2) > 6.188, 1), :);
%t_now = poi_local(end, 2)+0.00001;
%id = 13;
%poi_local = poi_local(poi_local(:,1)==id, :);
%v_no_reset = pm.ps * sum(exp(- 0.05*(t_now - poi_local(:, 2))))

l_min = min([length(ras_poi), length(ras_ref)]);
for k = 1 : l_min
  if any(ras_poi(k, :) - ras_ref(k, :))
    fprintf('first diff: poi=(%d, %.17g), ref=(%d, %.17g)\n',
      ras_poi(k, :), ras_ref(k, :));
    break;
  end
end
if k == l_min
  fprintf("Result: Correct. (Spike time exact match)\n");
end
ras_poi(1:5, :)
ras_ref(1:5, :)

return

fprintf('  ISI = %g\n', mean(ISI_poi));
fprintf('--> Result: Max diff poi (should = 0) = %g\n', maxabs(X_poi - X));


% vim: et sw=4 sts=4
