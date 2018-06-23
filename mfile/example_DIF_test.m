%

pm = [];
pm.neuron_model = 'DIF-single-GH';
pm.simu_method = 'simple';
pm.net  = zeros(1);
pm.nI   = floor(length(pm.net)/2);
pm.scee = 0.001;
pm.scie = 0.00;
pm.scei = 0.00;
pm.scii = 0.00;
pm.t    = 1e2;
pm.dt   = 1.0/32;
pm.stv  = pm.dt;
pm.seed = 123;
pm.pr   = 0;
pm.ps   = 0.01;
pm.pri  = 0;
pm.psi  = 0.012;
pm.extra_cmd = '';

% Suppose there are events to neuron 1 at time 5.0ms with strength 
% 0.01 (Excitatory) and 0.012 (negative means Inhibitory).
t_event_start = 5.0;
pm.input_event = [...
1 t_event_start pm.ps
1 t_event_start -pm.psi
];

% Set neuronal dynamical constants.
% For DIF-single-GH that is:
%   V_threshold, V_reset, V_leakage, V_excitatory, V_inhibitory, G_leak, synaptic_alpha,
%   tau_gE, tau_gE_s1, tau_gI, tau_gI_s1, Time_Refractory
% See single_neuron_dynamics.h and search "int n_var" for the details of other
%   neuron models.
% Currently, only LIF-GH, DIF-single-GH and HH-GH model support this function.
% And it can not combine with delayed synaptic or current/sine input.
pm.neuron_const = [...
1.0, 0.0, 0.0, 14/3.0, -2/3.0, 0.05, 1, 2.0, 0.5, 5.0, 0.8, 2.0
];

% Set initial state of the neurons.
% For DIF-single-GH that is
%   V, gE, gI, gE_s1, gI_s1
% The s1 means the equation that receives delta input directly.
% See single_neuron_dynamics.h and search "int n_var" for the details of other
%   neuron models.
% Note: no error checking is performed.
pm.initial_state = [...
0 0 0 0 0
];

[X, ISI, ras, pm1, extra_data] = gen_neu(pm, 'rm');

s_t = (1:length(X))*pm.stv;
figure(11);
plot(s_t, X);

f_exp_rise_fall = @(t, tg, th) tg*th*(exp(-t/tg) - exp(-t/th))/(tg-th) .* (t>0);
GE_ref = pm.ps  * f_exp_rise_fall(s_t-t_event_start, pm.neuron_const(8), pm.neuron_const(9));
GI_ref = pm.psi * f_exp_rise_fall(s_t-t_event_start, pm.neuron_const(10), pm.neuron_const(11));
G_ref = [GE_ref; GI_ref];

figure(21);
plot(s_t, extra_data.G, s_t, G_ref);

fprintf('Should be machine eps: \n');
maxabs(G_ref - extra_data.G)

