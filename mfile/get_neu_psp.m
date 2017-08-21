% Return an estimation of EPSP IPSP and the peak time.
% Usage:
%  % Ex1:
%  pm = [];
%  pm.prog_path = '../bin/gen_neu';
%  pm.neuron_model = 'HH-GH-cont-syn';
%  PSP = get_neu_psp(pm)
%
%  % Ex2:
%  PSP = get_neu_psp('HH-GH-cont-syn')

function PSP = get_neu_psp(pm0)
[~, tmp_f_name] = fileparts(tempname('./'));
events_file_path = ['./data/._tmp_neu_psp_poisson_' tmp_f_name '.txt'];

pm = [];
if ischar(pm0)
  pm.neuron_model = pm0;
else
  pm.neuron_model = pm0.neuron_model;
end
if isfield(pm0, 'prog_path')
  pm.prog_path = pm0.prog_path;
end
pm.simu_method  = 'auto';
pm.extra_cmd = sprintf('--input-event-path "%s"', events_file_path);

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

t_e = 50.0;          % poisson event time
pm.net = 'net_2_1';
pm.t  = 200 + t_e;   % assume the PSP will not last too long
pm.dt = 1/128.0;
pm.stv = pm.dt;
pm.pr = 0;
pm.ps = 1e-6;  % some very small value

fd = fopen(events_file_path, 'w');
fprintf(fd, '1 %.16e %.16e\n', t_e, pm.ps);
fprintf(fd, '2 %.16e %.16e\n', t_e, -pm.ps);
fclose(fd);

X = gen_neu(pm, 'new,rm', ['./data/.get_neu_psp_tmp' tmp_f_name]);
V_rest = X(1, floor(t_e/pm.stv) - 1);
X(:, 1:floor(t_e/pm.stv)) = [];
X = volt_unit * (X - V_rest);

% Get ps for 1 mV EPSP
% TODO: add interpolation to increase accuracy
[v_psp, pos_psp] = max(X(1,:));
PSP.mV_ps = pm.ps / v_psp;
PSP.t_ps = pos_psp*pm.stv;

% Get psi for 1 mV IPSP
[v_psp, pos_psp] = min(X(2,:));
PSP.mV_psi = -pm.ps / v_psp;
PSP.t_psi = pos_psp*pm.stv;

if v_psp == 0
  error('Fail to get EPSP (EPSP=0).');
end

switch pm.neuron_model
  case {'LIF-G', 'LIF-GH', 'HH-G', 'HH-GH'}
    PSP.mV_scee = PSP.mV_ps;
    PSP.t_scee  = PSP.t_ps;
    PSP.mV_scei = PSP.mV_psi;
    PSP.t_scei  = PSP.t_psi;
    delete(events_file_path);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the PSP for spike interaction

pm.net = 'net_2_2_T';
pm.t  = 300;
pm.dt = 1/128.0;
pm.stv = pm.dt;

pm.pr = 0;
pm.ps = 1.0 * PSP.mV_ps;

% Write events for neuron 2
fd = fopen(events_file_path, 'w');
t = t_e;
n_event = 0;
while (t < pm.t && n_event < 30)
  fprintf(fd, '2 %.16e\n', t);
  t = t + 0.3 + exp(-t*0.1);
  n_event = n_event + 1;
end
fclose(fd);

%% test scee
pm.nI = 0;
pm.scee = 1e-6;   % some small value
[X, ~, ras] = gen_neu(pm, 'new,rm', ['./data/.get_neu_psp_tmp' tmp_f_name]);
if isempty(ras)
  error('Fail to generate spike');
end
if size(ras,1) > 1
  warning('Too many spikes');
end
X(:, 1:floor(ras(1, 2)/pm.stv)) = [];
X = volt_unit * (X - V_rest);

[v_psp, pos_psp] = max(X(1, :));
PSP.mV_scee = pm.scee / v_psp;
PSP.t_scee = pos_psp*pm.stv;

%% test scei
pm.nI = 1;
pm.scei = 1e-6;
[X, ~, ras] = gen_neu(pm, 'new,rm', ['./data/.get_neu_psp_tmp' tmp_f_name]);
X(:, 1:floor(ras(1, 2)/pm.stv)) = [];
X = volt_unit * (X - V_rest);

[v_psp, pos_psp] = min(X(1, :));
PSP.mV_scei = -pm.scei / v_psp;
PSP.t_scei = pos_psp*pm.stv;

delete(events_file_path);

end
