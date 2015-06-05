%

clear('pm');          % a new parameter set, or pm = [];
pm.prog_path = '../bin/gen_neu';
pm.neuron_model = 'LIF-G';  % one of LIF-G, LIF-GH, HH-GH
pm.net  = 'net_1_0';  % can also be a connectivity matrix or full file path
pm.nI   = 0;          % default: 0. Number of Inhibitory neurons.
                      %             Indexes are later half
pm.scee = 0.006;
pm.scie = 0.00;       % default: 0. Strength from Ex. to In.
pm.scei = 0.00;       % default: 0. Strength from In. to Ex.
pm.scii = 0.00;       % default: 0.
pm.pr   = 1.6;        % can be a vector (not yet)
pm.ps   = 0.006;       % can be a vector (not yet)
pm.t    = 1e5;
pm.dt   = 1.0/32;     % default: 1/32
pm.stv  = 0.5;        % default: 0.5
pm.seed = 'auto';     % default: 'auto'(or []). Accept integers
pm.extra_cmd = '';    % default: '--RC-filter 0 1'
                      % put all other parameters here.
%[X, ISI, ras] = gen_neu(pm, 'cmd');
tic;
[X, ISI, ras] = gen_neu(pm, 'new');
toc;

ISI
figure(1);
hist(diff(ras(ras(:,1)==1, 2)), 200);
xlim([0, 60]);

clear('pm');          % a new parameter set, or pm = [];
pm.prog_path = '../bin/gen_neu';
pm.neuron_model = 'LIF-GH';  % one of LIF-G, LIF-GH, HH-GH
pm.net  = 'net_1_0';  % can also be a connectivity matrix or full file path
pm.nI   = 0;          % default: 0. Number of Inhibitory neurons.
                      %             Indexes are later half
pm.scee = 0.013;
pm.scie = 0.00;       % default: 0. Strength from Ex. to In.
pm.scei = 0.00;       % default: 0. Strength from In. to Ex.
pm.scii = 0.00;       % default: 0.
pm.pr   = 1.6;        % can be a vector (not yet)
pm.ps   = 0.012;       % can be a vector (not yet)
pm.t    = 1e5;
pm.dt   = 1.0/32;     % default: 1/32
pm.stv  = 0.5;        % default: 0.5
pm.seed = 'auto';     % default: 'auto'(or []). Accept integers
pm.extra_cmd = '';    % default: '--RC-filter 0 1'
                      % put all other parameters here.
%[X, ISI, ras] = gen_neu(pm, 'cmd');
tic;
[X, ISI, ras] = gen_neu(pm, 'new');
toc;

ISI
figure(2);
hist(diff(ras(ras(:,1)==1, 2)), 200);
xlim([0, 60]);

clear('pm');          % a new parameter set, or pm = [];
pm.prog_path = '../bin/gen_neu';
pm.neuron_model = 'HH-GH';  % one of LIF-G, LIF-GH, HH-GH
pm.net  = 'net_1_0';  % can also be a connectivity matrix or full file path
pm.nI   = 0;          % default: 0. Number of Inhibitory neurons.
                      %             Indexes are later half
pm.scee = 0.05;
pm.scie = 0.00;       % default: 0. Strength from Ex. to In.
pm.scei = 0.00;       % default: 0. Strength from In. to Ex.
pm.scii = 0.00;       % default: 0.
pm.pr   = 1.6;        % can be a vector (not yet)
pm.ps   = 0.04;       % can be a vector (not yet)
pm.t    = 1e5;
pm.dt   = 1.0/32;     % default: 1/32
pm.stv  = 0.5;        % default: 0.5
pm.seed = 'auto';     % default: 'auto'(or []). Accept integers
pm.extra_cmd = '';    % default: '--RC-filter 0 1'
                      % put all other parameters here.
%[X, ISI, ras] = gen_neu(pm, 'cmd');
tic;
[X, ISI, ras] = gen_neu(pm, 'new');
toc;

ISI
figure(3);
hist(diff(ras(ras(:,1)==1, 2)), 200);
xlim([0, 60]);
%plot(X');
