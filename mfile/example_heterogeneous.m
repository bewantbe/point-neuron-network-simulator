%

pm = [];
pm.neuron_model = 'LIF-GH';
pm.simu_method = 'simple';
pm.net     = zeros(1);
pm.nI      = floor(length(pm.net)/2);
pm.scee_mV = 0.001;
pm.scie_mV = 0.00;
pm.scei_mV = 0.00;
pm.scii_mV = 0.00;
pm.t    = 1e2;
pm.dt   = 1.0/32;
pm.stv  = pm.dt;
pm.seed = 123;
pm.pr    = 0;
pm.ps    = 0.01;
pm.pri   = 0;
pm.psi   = 0.012;
pm.extra_cmd = '';
pm.input_event = [...
1 5.0 0.01
1 5.0 -0.012
];
pm.tau_g = [...
1.0 2.0 1.0 2.0
];
% Set initial state of the neurons.
% For LIF-GH that is V, gE, gI, gE_s1, gI_s1
%   The s1 means the equation that receives delta input directly.
% See single_neuron_dynamics.h and search "int n_var" for the details of other
%   neuron models.
% Note: no error checking is performed.
pm.initial_state = [...
0 0 0 0 0
];

[X, ISI, ras, ~, extra_data] = gen_neu(pm, 'rm');

f_exp_rise_fall = @(t, tg, th) tg*th*(exp(-t/tg) - exp(-t/th))/(tg-th) .* (t>0);

s_t = (1:length(X))*pm.stv;

GE_ref = pm.ps  * f_exp_rise_fall(s_t-5.0, pm.tau_g(1), pm.tau_g(2));
GI_ref = pm.psi * f_exp_rise_fall(s_t-5.0, pm.tau_g(3), pm.tau_g(4));
G_ref = [GE_ref; GI_ref];

figure(20);
plot(s_t, extra_data.G, s_t, G_ref);

fprintf('Should be machine eps: \n');
maxabs(G_ref - extra_data.G)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Long time test.
pm = [];
pm.neuron_model = 'LIF-GH';
pm.simu_method = 'simple';
pm.net     = ones(100);
pm.nI      = floor(length(pm.net)/2);
pm.scee_mV = 0.001;
pm.scie_mV = 0.00;
pm.scei_mV = 0.00;
pm.scii_mV = 0.00;
pm.t    = 1e3;
pm.dt   = 1.0/32;
pm.stv  = pm.dt;
pm.seed = 123;
pm.pr    = 1;
pm.ps    = 0.01;
pm.pri   = 1;
pm.psi   = 0.012;
pm.extra_cmd = '';

[X_ref, ISI_ref, ras_ref] = gen_neu(pm, 'rm');

pm.tau_g = ones(size(pm.net,1), 1) * [0.5 2.0 0.8 5.0];
[X, ISI, ras] = gen_neu(pm, 'rm');

fprintf('Should be exactly "0":\n');
maxabs(X-X_ref)

