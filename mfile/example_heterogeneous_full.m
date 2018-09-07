% Demo: How to set every dynamical variables in the neuronal model.

pm = [];
pm.neuron_model = 'LIF-GH';
pm.simu_method = 'big-delay';
pm.net  = 1*(eye(2)==0);
pm.nI   = floor(length(pm.net)/2);
pm.scee = 0.03;
pm.scie = 0.03;
pm.scei = 0.05;
pm.scii = 0.05;
pm.t    = 1e2;
pm.dt   = 1.0/32;
pm.stv  = pm.dt;
pm.seed = 123;
pm.pr   = 0;
pm.ps   = 0.02;
pm.pri  = 0;
pm.psi  = 0.03;
pm.extra_cmd = '';
pm.synaptic_delay = 5;

% Suppose there are events to neuron 1 at time 5.0ms with strength 
% 0.01 (Excitatory) and 0.012 (negative means Inhibitory).
pm.input_event = [...
1 5.0 0.01
1 5.0 -0.012
];

% Set neuronal dynamical constants.
% For LIF-GH that is:
%   V_threshold, V_reset, V_leakage, V_excitatory, V_inhibitory, G_leak, tau_gE,
%   tau_gE_s1, tau_gI, tau_gI_s1, Time_Refractory
% See single_neuron_dynamics.h and search "int n_var" for the details of other
%   neuron models.
% Currently, only LIF-GH and HH-GH model support this function. And it can not 
%   combine with heterogeneous delayed synaptic or current/sine input.
pm.neuron_const = ones(length(pm.net), 1) * [...
1.0, 0.0, 0.0, 14/3.0, -2/3.0, 0.05, 2.0, 0.5, 5.0, 0.8, 2.0
];

% Set initial state of the neurons.
% For LIF-GH that is V, gE, gI, gE_s1, gI_s1
%   The s1 means the equation that receives delta input directly.
% See single_neuron_dynamics.h and search "int n_var" for the details of other
%   neuron models.
% Note: no error checking is performed.
pm.initial_state = ones(length(pm.net), 1) * [...
0 0 0 0 0
];

[X, ISI, ras, pm1, extra_data] = gen_neu(pm, 'rm');

s_t = (1:length(X))*pm.stv;
figure(10);
plot(s_t, X);

f_exp_rise_fall = @(t, tg, th) tg*th*(exp(-t/tg) - exp(-t/th))/(tg-th) .* (t>0);
GE_ref = pm.input_event(1,3) * f_exp_rise_fall(s_t-pm.input_event(1,2), pm.neuron_const(1,7), pm.neuron_const(1,8));
GI_ref =-pm.input_event(2,3) * f_exp_rise_fall(s_t-pm.input_event(2,2), pm.neuron_const(1,9), pm.neuron_const(1,10));
G_ref = [GE_ref; GI_ref];

figure(20);
plot(s_t, extra_data.G(1:2, :), s_t, G_ref);

fprintf('Should be machine eps: \n');
maxabs(G_ref - extra_data.G(1:2, :))

%% Experiment 2, test the interaction.
pm.input_event = zeros(25, 2);
pm.input_event(:,1) = 1;
pm.input_event(:,2) = 1:25;  % stimulus at 1, 2, ..., 25 ms
[X, ISI, ras] = gen_neu(pm, 'rm');

figure(110);
plot(s_t, X);

id_up = find(X(2,:)>0, 1);
pp = polyfit(0:2, X(2,id_up:id_up+2),2);
t_rec = fzero(@(x) polyval(pp, x), [-1,0]);
t_rec = (t_rec + id_up) * pm.stv;
fprintf('synaptic delay = %g ms\n', t_rec - ras(1,2));

