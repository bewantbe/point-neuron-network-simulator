% Test consistency across versions

function search_defe(k_case)

addpath('../mfile');
maxabs = @(x) max(abs(x(:)));

path_tag_executable = '../bin/gen_neu';
session_random = char(randi(26, 1, 10) - 1 + 'a');

for seed = (1:2000)+2000*k_case
    pm = [];
    pm.prog_path = path_tag_executable;
    pm.neuron_model = 'LIF-G';
    pm.simu_method  = 'auto';
%    pm.synaptic_delay = 0.5;
%    pm.simu_method  = 'SSC';
    pm.net     = ones(2);
    pm.nI      = 0;
    pm.scee_mV = 1.0;
    pm.pr      = 4.0;
    pm.ps_mV   = 1.0;
    pm.t    = 1e2;
    pm.dt   = 1.0/32;
    pm.stv  = 0.5;
    pm.seed = seed;
    pm.extra_cmd = '';

    [~, tmp_f_name] = fileparts(tempname('./'));
    poi_file = ['./data/.poi_' tmp_f_name '.txt'];

    % Test new version, reproduce, output of poisson events
    pm.extra_cmd = ['--output-poisson-events-path ' poi_file];
    [X, ISI, ras, pm] = gen_neu(pm, 'new,rm', ['./data/' session_random]);

    % construct input events to neuron "id_test".
    id_test = 1;
    % input from poisson
    poi = load(poi_file);
    poi = poi(poi(:,1)==id_test, :);  % keep only id_test related spikes
    poi(:,1) = 1;
    % input from other neurons
    ras = ras(ras(:,1)~=id_test, :);
    sc_strength = [pm.scee, pm.scie, pm.scei, pm.scii];
    ras_strength = pm.net_adj(id_test, ras(:, 1)) .* sc_strength((id_test>pm.nE) + 2*(ras(:, 1)>pm.nE) + 1);    
    ras = [ras ras_strength'];
    ras(ras(:, 3)==0, :) = [];
    ras(:, 1) = 1;
    
%        if isfield(pm, 'synaptic_delay')
%            ras(:, 2) += pm.synaptic_delay;
%        end
%        
    poi = [poi; ras];
    [~, id_sort] = sort(poi(:,2));
    poi = poi(id_sort, :);
    
    pm.input_event = poi;
    pm.extra_cmd = '';
    pm.net     = zeros(1);
    pm.nE      = 1;
    pm.nI      = 0;
    X_single = gen_neu(pm, 'new,rm', ['./data/' session_random]);
    
    delete(poi_file);

    % The error is generally not zero, but will converge to zero when
    % dt -> 0. The problem is that the Spike-Correction algorithm (SSC)
    % sub-divides the time intervals more than necessary.
%    fprintf('--> Result: Max diff single = %g\n', maxabs(X(1,:) - X_single));

%    d = X(1,1:200) - X_single(1:200);
    d = X(1,:) - X_single;
    if maxabs(d) > 0
        fprintf('%d --> Max diff = %g @%d ', seed, maxabs(d), find(d, 1));
    else
        fprintf('%d --> Max diff = %g ', seed, maxabs(d));
    end
    fprintf(' ISI=%.16g\n', mean(ISI));
end

% vim: et sw=4 sts=4
