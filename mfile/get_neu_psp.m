% Return an estimation of EPSP IPSP and the peak time.
% Usage:
%  pm = [];
%  pm.prog_path = '../bin/gen_neu';
%  pm.neuron_model = 'HH-GH-cont-syn';
%  PSP = get_neu_psp(pm)

function PSP = get_neu_psp(pm)
events_file_path = 'poisson_events.txt';

pm.simu_method  = 'auto';
pm.extra_cmd = sprintf('--input-event-path %s', events_file_path);

volt_unit = 1.0;
switch pm.neuron_model
  case 'HH-GH-cont-syn'
    volt_unit = 10.0;
  case {'LIF-G' 'LIF-GH'}
    volt_unit = 15.0;
end

PSP.neuron_model = pm.neuron_model;
PSP.volt_unit = volt_unit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the PSP for external input

t_e = 50.0;

fd = fopen(events_file_path, 'w');
fprintf(fd, '0 %.16e\n', t_e);
fclose(fd);

pm.net = 'net_1_0';
pm.t  = 200 + t_e;   % assume the PSP will not last too long
pm.dt = 1/128.0;
pm.stv = pm.dt;

pm.pr = 0;
pm.ps = 1e-6;  % some very small value

X = gen_neu(pm, 'new,rm');
V_rest = X(1, floor(t_e/pm.stv) - 1);
X(:, 1:floor(t_e/pm.stv)) = [];
X = volt_unit * (X - V_rest);

[v_psp, pos_psp] = max(X);
t_psp = pos_psp*pm.stv;

mV_EPSP_ps = pm.ps / v_psp;

PSP.mV_ps = mV_EPSP_ps;
PSP.t_ps = t_psp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the PSP for spike interaction

pm.net = 'net_2_2_T';
pm.t  = 500;
pm.dt = 1/128.0;
pm.stv = pm.dt;

pm.pr = 0;
pm.ps = 1.0 * mV_EPSP_ps;

% write events for neuron 2
fd = fopen(events_file_path, 'w');
t = t_e;
n_event = 0;
while (t < pm.t && n_event < 30)
  fprintf(fd, '1 %.16e\n', t);
  t = t + 0.3 + exp(-t*0.1);
  n_event = n_event + 1;
end
fclose(fd);

%% test scee

pm.nI = 0;
pm.scee = 1e-6;   % some small value

[X, ~, ras] = gen_neu(pm, 'new,rm');
if isempty(ras)
  error('failed to generate spike');
end
if size(ras,1) > 1
  error('too many spikes');
end

X(:, 1:floor(ras(1, 2)/pm.stv)) = [];
X = volt_unit * (X - V_rest);

[v_psp, pos_psp] = max(X(1, :));
t_psp = pos_psp*pm.stv;

mV_EPSP_scee = pm.scee / v_psp;
PSP.mV_scee = mV_EPSP_scee;
PSP.t_scee = t_psp;

%% test scei

pm.nI = 1;
pm.scei = 1e-6;
[X, ~, ras] = gen_neu(pm, 'new,rm');
X(:, 1:floor(ras(1, 2)/pm.stv)) = [];
X = volt_unit * (X - V_rest);
[v_psp, pos_psp] = min(X(1, :));
t_psp = pos_psp*pm.stv;

mV_IPSP_scei = -pm.scei / v_psp;
PSP.mV_scei = mV_IPSP_scei;
PSP.t_scei = t_psp;

delete(events_file_path);

end
