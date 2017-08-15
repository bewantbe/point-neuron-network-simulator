% Test consistency across versions
addpath('../mfile');
maxabs = @(x) max(abs(x(:)));

s_neuron_model = {'LIF-G', 'LIF-GH', 'HH-G', 'HH-GH', 'HH-GH-cont-syn'};
path_ref_executable = '../bin/gen_neu_ref';
path_tag_executable = '../bin/gen_neu';

id_nm = 1;

for id_nm = 1#:length(s_neuron_model)
    fprintf('===== Testing model: %s =====\n', s_neuron_model{id_nm});
    pm = [];
    pm.prog_path = path_ref_executable;
    pm.neuron_model = s_neuron_model{id_nm};
    pm.simu_method  = 'auto';
    pm.net     = ones(15);
    pm.nI      = 0;
    pm.scee_mV = 1.0;
    pm.pr      = 1.0;
    pm.ps_mV   = 0.9;
    pm.t    = 1e4;
    pm.dt   = 1.0/32;
    pm.stv  = 0.5;
    pm.seed = 324789;
    pm.extra_cmd = '-v';

    [X_ref, ISI] = gen_neu(pm, 'new,rm');
    
    pm.prog_path = path_tag_executable;
    pm.extra_cmd = [pm.extra_cmd ' --output-poisson-events-path poi.txt'];
    [X, ~, ras, pm] = gen_neu(pm, 'new,rm');

    % construct input events to neuron "id_test".
    id_test = 1;
    % input from poisson
    poi = load('poi.txt');
    poi = poi(poi(:,1)==id_test, :);  % keep only id_test related spikes
    poi(:,1) = 1;
    % input from other neurons
    ras = ras(ras(:,1)~=id_test, :);
    sc_strength = [pm.scee, pm.scie, pm.scei, pm.scii];
    ras_strength = pm.net_adj(id_test, ras(:, 1)) .* sc_strength((id_test>pm.nE) + 2*(ras(:, 1)>pm.nE) + 1);    
    ras = [ras ras_strength'];
    ras(ras(:, 3)==0) = [];
    ras(:, 1) = 1;
    
    poi = [poi; ras];
    [~, id_sort] = sort(poi(:,2));
    poi = poi(id_sort, :);
    
    pm.input_event = poi;
    pm.extra_cmd = '-v';
    pm.net     = zeros(1);
    pm.nE      = 1;
    pm.nI      = 0;
    X_single = gen_neu(pm, 'new,rm');

    fprintf('  ISI = %g\n', mean(ISI));
    fprintf('*** Result: Max diff = %g\n', maxabs(X - X_ref));
    fprintf('*** Result: Max diff single = %g\n', maxabs(X(1,:) - X_single));
end

%plot(X, X_ref)
rg = 1:200;
figure(1);
plot(X(1, rg) - X_single(1, rg))

figure(2);
plot(X(1, rg))

% vim: et sw=4 sts=4
