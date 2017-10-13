% Neuron network simulator (interface to gen_neu).
% Can use cached data automatically.
%
%  [X, isi, ras, pm, extra_data] = gen_neu(pm [, gen_cmd [, data_dir_prefix]])
%
% Usage example 1:       % the items with default value are optional
%  pm = [];
%  %pm.prog_path = '../bin/gen_neu';  % path to the executable.
%  pm.neuron_model = 'HH-GH';  % LIF-G, LIF-GH, HH-GH etc. See gen_neu --help
%  pm.simu_method = 'SSC';     % Spike-Spike-Correction, 'simple' or 'auto'
%  pm.net  = [0 1; 0 0]; % Connectivity matrix or file to the matrix.
%  pm.nI   = 0;          % default: 0. Number of Inhibitory neurons.
%                        %             Indexes are later half.
%  pm.scee_mV = 0.5;
%  pm.scie_mV = 0.0;     % default: 0. Strength from Ex. to In.
%  pm.scei_mV = 0.0;     % default: 0. Strength from In. to Ex.
%  pm.scii_mV = 0.0;     % default: 0.
%  pm.pr      = 1.6;     % poisson rate, can be a vector
%  pm.ps_mV   = 0.4;     % poisson strength, can be a vector
%  pm.t    = 1e4;
%  pm.dt   = 2^-5;       % default: 1/32
%  pm.stv  = 0.5;        % default: 0.5
%  pm.seed = 'auto';     % default: 'auto'. Accept an integer or an integer vector
%  pm.extra_cmd = '-v';  % All other parameters go here.
%  [X, ISI, ras] = gen_neu(pm);
%
% Usage example 2: Always re-generate data, then read it.
%  X = gen_neu(pm, 'new');
%
% Usage example 3: Generate data if not exist, read it, then remove data files.
%  X = gen_neu(pm, 'rm');
%
% Usage example 4: Get also conductance and gating variables.
%  [X, isi, ras, pm, extra_data] = gen_neu(pm,'new,extra_data');
%  extra_data.G(1)  % E conductance of #1 neuron.
%  extra_data.G(2)  % I conductance of #1 neuron.
%  extra_data.G(3)  % E conductance of #2 neuron.
%  extra_data.G(4)  % I conductance of #2 neuron.
%  etc.
%
% Other possible values for "gen_cmd":
%  'legancy'  Call the old good raster_tuning.
%  'read'   Read data files if exist, otherwise do nothing and return [];
%  'nameX'  return path to the voltage data file, instead of read it;
%  'cmd'    Show command call to raster_tuning_HH, then exit. Useful for debug;
%  'ext_T'  Generate a bit more data, so reduce "head effect";
%  'v'      Be verbose;
%  'h'      Show help for the command line generator.
%
% Return value pm is a "normalized" input parameter set.

% Behaviour (parameter gen_cmd)
% gen_cmd   no data    have data   background(partial data)
% ''        gen+read   read        gen+read(overwrite)
% 'new'     gen+read   gen+read    gen+read(overwrite)
% 'nameX'   gen+name   name        gen+read(overwrite)
% 'read'    none       read        exit
% 'rm'      none       rm          none
% 'cmd'     print cmd for gen then exit
% 'ext_T'
% 'extra_data'
% background  gen&     none        gen&

function [X, ISI, ras, pm, extra_data] ...
         = gen_neu(pm, gen_cmd, data_dir_prefix)
t_start = tic;
if nargin()==0
    disp(' [X, ISI, ras] = gen_neu(pm [, gen_cmd [, data_dir_prefix]])');
    disp(' Type "help gen_neu" for more help');
    error('No input parameter.');
end
if ~exist('fflush', 'builtin')
    fflush = @(a) 0;
    stdout = 0;
end
X=[];
ISI=[];
ras=[];

if ~ has_nonempty_field(pm, 'prog_path')
    % Search gen_neu in "../bin" then "./"
    pathdir = fileparts(mfilename('fullpath'));
    exepath = sprintf('%s%s..%sbin%sgen_neu', pathdir, filesep, filesep, filesep);
    if ~exist(exepath, 'file') && ~exist([exepath '.exe'], 'file')
        exepath = sprintf('%s%sgen_neu', pathdir, filesep);
    end
    if ~exist(exepath, 'file') && ~exist([exepath '.exe'], 'file')
        error('Executable gen_neu not found!');
    end
    pm.prog_path = exepath;
end

if ~isfield(pm, 'extra_cmd')
    pm.extra_cmd = '';
end

% Default generator settings
new_run        = false;
return_X_name  = false;
mode_rm_only   = false;
mode_show_cmd  = false;
mode_read_only = false;
mode_legancy   = ~isempty(strfind(pm.prog_path, 'raster_tuning'));
mode_extra_data = nargout >= 5;
mode_run_in_background = ~isempty(find(pm.extra_cmd == '&', 1, 'last'));

ext_tmp = [' ' strtrim(pm.extra_cmd) ' '];
b_verbose = ~isempty([strfind(ext_tmp,' -v ') strfind(ext_tmp,' --verbose ')]);  % show time cost
ext_T = 0;

% Read generator parameters
if ~exist('gen_cmd','var') || isempty(gen_cmd)
    gen_cmd = '';
end
gen_cmd = strtrim(gen_cmd);
while ~isempty(gen_cmd)
    [tok, gen_cmd] = strtok(gen_cmd, ' ,');
    switch tok
    case 'rm'
        mode_rm_only = true;    % Remove specific data files then exit
    case 'read'
        mode_read_only = true;  % Read data files if there are, otherwise exit
    case 'new'
        new_run = true;         % Regenerate data and read it
    case 'nameX'
        return_X_name = true;   % Return file path of voltage data file
    case 'legancy'
        mode_legancy  = true;   % Show the command to call, then exit
    case 'cmd'
        mode_show_cmd = true;   % Show the command to call, then exit
    case 'ext_T'
        ext_T = 1e2;            % extra calculation time
    case 'extra_data'
        mode_extra_data = true; % output extra data: G, h, m, n
        % data will be in struct extra_data
    case {'verbose', 'v'}
        b_verbose = true;
    case 'h'
        system(['"' pm.prog_path '" -h']);
        return
    otherwise
        error('no this option: "%s"', tok);
    end
end

if ~exist('data_dir_prefix', 'var')
    data_dir_prefix = ['.', filesep, 'data', filesep];
end

% Fields related to strength.
field_v = {'scee', 'scie', 'scei', 'scii', 'ps', 'psi'};

% Test mis-spelling.
for fv = field_v
    if isfield(pm, [fv{1} '_mv'])
        error(['mis-spelling: ' fv{1} '_mv should be ' fv{1} '_mV']);
    end
end

% Convert "mV" to internal unit.
% if contains field that needs conversion
if any(cellfun(@(fv) isfield(pm, [fv '_mV']), field_v))
    PSP = get_neu_psp(pm);
    PSP_v = [PSP.mV_scee, PSP.mV_scee, PSP.mV_scei, PSP.mV_scei, PSP.mV_ps, PSP.mV_psi];
    for id_fv = 1:length(field_v)
        fv = field_v{id_fv};
        if has_nonempty_field(pm, [fv '_mV'])
            s = pm.([fv '_mV']) * PSP_v(id_fv);
            if isfield(pm, fv) && s ~= pm.(fv)
                error(['incompatible EPSP/IPSP strength specification: ' fv ' ' fv '_mV']);
            end
            pm.(fv) = s;
        end
    end
end

% Fill-in default values.
if ~mode_legancy
    if ~ has_nonempty_field(pm, 'neuron_model')
        warning('gen_neu:model', 'pm.neuron_model not specified! Using "HH-PT-GH"');
        pm.neuron_model = 'HH-PT-GH';
    end
    if ~ has_nonempty_field(pm, 'simu_method')
        warning('gen_neu:model', 'pm.simu_method not specified! Using "auto".');
        pm.simu_method = 'auto';
    end
else
    if ~ has_nonempty_field(pm, 'neuron_model')
        [~, pm.neuron_model] = fileparts(pm.prog_path);
        % Try to convert program name to known model.
        if strmatch('raster_tuning_HH', pm.neuron_model)
            pm.neuron_model = 'legancy-HH-GH-cont-syn';
        elseif strmatch('raster_tuning_LIF_GH', pm.neuron_model)
            pm.neuron_model = 'legancy-LIF-GH';
        elseif strmatch('raster_tuning_LIF', pm.neuron_model)
            pm.neuron_model = 'legancy-LIF-G';
        end
    end
    if ~b_verbose
        pm.extra_cmd = ['-q' pm.extra_cmd];
    end
end
if ~ has_nonempty_field(pm, 'stv')
    pm.stv = 0.5;
end
if ~ has_nonempty_field(pm, 'dt')
    pm.dt = 1.0/32;
end
if ~ has_nonempty_field(pm, 'scee')
    pm.scee = 0;
end
if ~ has_nonempty_field(pm, 'scie')
    pm.scie = 0;  % Strength from Ex. to In.
end
if ~ has_nonempty_field(pm, 'scei')
    pm.scei = 0;  % Strength from In. to Ex.
end
if ~ has_nonempty_field(pm, 'scii')
    pm.scii = 0;
end
if ~ has_nonempty_field(pm, 'seed')
    if ~mode_legancy
        pm.seed = randi(2^32-1, 1, 2);
    else
        pm.seed = randi(2^32-1, 1, 1);
    end
end
if has_nonempty_field(pm, 't_warming_up')
    ext_T = 0;  % use pm.t_warming_up instead of ext_T
end
if ext_T > 0
    pm.ext_T = ext_T;
end

if xor(isfield(pm,'pr'), isfield(pm,'ps'))
    warning('gen_neu:pm', 'pr and ps not privided together.');
end
if xor(isfield(pm,'pri'), isfield(pm,'psi'))
    warning('gen_neu:pm', 'pri and psi not privided together.');
end

if ~ has_nonempty_field(pm, 'net')
    pm.net = 'net_1_0';
end

if ~ischar(pm.net)
  b_clean_net_file = true;
else
  b_clean_net_file = false;
end

% Prepare the neuronal network.
if ischar(pm.net)
    [network, mat_path] = get_network(pm.net, data_dir_prefix);
    if isempty(network)
        if has_nonempty_field(pm, 'net_adj') && isnumeric(pm.net_adj) ...
           && diff(size(pm.net_adj)) == 0
            network = pm.net_adj;
            [mat_path, pm.net] = save_network(pm.net_adj, data_dir_prefix);
            fprintf('Using pm.net_adj as the network\n');
        else
            error('Set pm.net to the file name of net or adjacency matrix');
        end
    end
    [~, pm.net] = fileparts(pm.net);          % Use the name without extension
else
    % so pm.net is connectivity matrix?
    % save this matrix, so that it can be read by `gen_neu'
    network = pm.net;
    [mat_path, pm.net] = save_network(pm.net, data_dir_prefix);
end
pm.net_path = mat_path;
pm.net_adj  = network;

p = size(network,1);
if ~ has_nonempty_field(pm, 'nI')
    pm.nI = 0;  % number of inhibitory neurons
end
if ~ has_nonempty_field(pm, 'nE')
    pm.nE = p - pm.nI;
end
if pm.nI + pm.nE ~= p || pm.nE < 0 || pm.nI < 0
    fprintf('  pm.nE(%d) + pm.nI(%d) = %d, p_net = %d\n', ...
        pm.nE, pm.nI, pm.nE + pm.nI, p);
    error('gen_neu: Number of neurons inconsist with the network!');
end

if has_nonempty_field(pm, 'synaptic_net_delay')
    net_delay_path = save_network(pm.synaptic_net_delay, ...
                     [data_dir_prefix 'net_delay_']);
else
    net_delay_path = '';
end

if has_nonempty_field(pm, 'input_event')
    [~, tmp_f_name] = fileparts(tempname('./'));
    poisson_path = [data_dir_prefix 'external_event_' tmp_f_name '.txt'];
    if (size(pm.input_event,2) ~= 2 && size(pm.input_event,2) ~= 3)
        error('pm.input_event must be 2 or 3 column, the order is "id time [strength]"');
    end
    func_fraction = @(x) x - floor(x);
    if any(func_fraction(pm.input_event(:,1)))
        error('pm.input_event first column must be neuron index(1 based).');
    end
    if ~issorted(pm.input_event(:,2))
        warning('gen_neu:pm', 'pm.input_event is not time sorted, sorting for you.');
        [~, id_sort] = sort(pm.input_event(:,2));
        pm.input_event = pm.input_event(id_sort, :);
    end
    fd = fopen(poisson_path, 'w');
    if size(pm.input_event, 2) == 2
        fprintf(fd, '%d %.16e\n', pm.input_event');
    else
        fprintf(fd, '%d %.16e %.16e\n', pm.input_event');
    end
    fclose(fd);
else
    poisson_path = '';
end

if has_nonempty_field(pm, 'force_spikes')
    [~, tmp_f_name] = fileparts(tempname(['.' filesep]));
    force_spike_path = [data_dir_prefix tmp_f_name 'force_spike.txt'];
    fd = fopen(force_spike_path, 'w');
    fprintf(fd, '%d %.17g\n', pm.force_spikes');
    fclose(fd);
else
    force_spike_path = '';
end

% construct file paths
pr_name = get_prps_str(pm, 'pr');
ps_name = get_prps_str(pm, 'ps');
pri_name = get_prps_str(pm, 'pri');
psi_name = get_prps_str(pm, 'psi');
st_sc = strrep(mat2str([pm.scee, pm.scie, pm.scei, pm.scii], 5),' ',',');
st_p  = strrep(mat2str([pm.nE, pm.nI]),' ',',');
file_inf_st =...
    sprintf('%s_p=%s_sc=%s_%s_%s_%s_%s_stv=%g_t=%.2e',...
            pm.net, st_p(2:end-1), st_sc(2:end-1), pr_name,...
            ps_name, pri_name, psi_name, pm.stv, pm.t);
file_prefix = [data_dir_prefix, pm.neuron_model, '_'];
output_name     = [file_prefix, 'volt_',file_inf_st,'.dat'];
output_ISI_name = [file_prefix, 'ISI_', file_inf_st,'.txt'];
output_RAS_name = [file_prefix, 'RAS_', file_inf_st,'.txt'];
output_G_name   = [file_prefix, 'G_',file_inf_st,'.dat'];
output_gating_name = [file_prefix, 'gating_',file_inf_st,'.dat'];

st_paths =...
    sprintf('--volt-path "%s" --isi-path "%s" --ras-path "%s"',...
            output_name, output_ISI_name, output_RAS_name);
if mode_extra_data
    st_paths = [st_paths ...
        sprintf(' --conductance-path "%s" --ion-gate-path "%s"',...
                output_G_name, output_gating_name)];
end

% Path string escape
f_unescape = @(s) strrep(strrep(s, '%', '%%'), '\', '\\');

if ~mode_legancy
    % Parameters that will pass to executable gen_neu.
    c_options = {...
    {'prog_path', [], ''}
    {'neuron_model', '%s'}
    {'simu_method', '%s', '--simulation-method'}
    {'nE'}
    {'nI'}
    {'net_path', [], '--net'}
    {'net_adj', '', inlineif(issparse(pm.net_adj), '--sparse-net', '')}
    {'pr'}
    {'ps'}
    {'pri'}
    {'psi'}
    {'scee'}
    {'scie'}
    {'scei'}
    {'scii'}
    {'extI'}
    {'sine_amp', [], '--current-sine-amp'}
    {'sine_freq', [], '--current-sine-freq'}
    {'t', sprintf('%.17g', pm.t + ext_T)}
    {'dt'}
    {'stv'}
    {'t_warming_up'}
    {'seed'}
    {'spike_threshold', [], '--set-threshold'}
    {'synaptic_delay'}
    {'synaptic_net_delay', '', ['--synaptic-net-delay "' net_delay_path '"']}
    {'input_event', '', ['--input-event-path "' poisson_path '"']}
    {'force_spikes', '', ['--force-spike-list "' force_spike_path '"']}
    {'prog_path', '', st_paths}
    {'output_poisson_path'}
    {'extra_cmd', '%s', ''}
    % Recognized parameters that will not pass to gen_neu.
    {'net', '', ''}
    {'cmd_str', '', ''}
    {'ext_T', '', ''}
    {'scee_mV', '', ''}
    {'scie_mV', '', ''}
    {'scei_mV', '', ''}
    {'scii_mV', '', ''}
    {'ps_mV', '', ''}
    {'psi_mV', '', ''}
    }.';
else
    if mode_extra_data
        pm.save_conductance = output_G_name;
    end
    if isfield(pm,'pr') && length(pm.pr) > 1
        pm.pr_mul = pm.pr;
        pm.pr = 1;
    end
    if isfield(pm,'ps') && length(pm.ps) > 1
        pm.ps_mul = pm.ps;
        pm.ps = 1;
    end
    if isfield(pm,'psi') && length(pm.psi) > 1
        pm.psi_mul = pm.psi;
        pm.psi = 1;
    end
    % Parameters that will pass to executable gen_neu.
    c_options = {...
    {'prog_path', '"%s" -ng -inf - --RC-filter 0 1', ''}
    {'nE', [], '-n'}
    {'nI', '%d', ''}
    {'net_path', [], '-mat'}
    {'pr', [], '-pr'}
    {'ps', [], '-ps'}
    {'psi', [], '-psi'}
    {'pr_mul'}
    {'ps_mul'}
    {'psi_mul'}
    {'scee', [], '-scee'}
    {'scie', [], '-scie'}
    {'scei', [], '-scei'}
    {'scii', [], '-scii'}
    {'t', sprintf('%.17g', pm.t + ext_T), '-t'}
    {'dt', [], '-dt'}
    {'stv', [], '-stv'}
    {'seed', [], '-seed'}
    {'prog_path', '', ['--bin-save -o "' output_name '"']}
    {'prog_path', '', ['--save-spike-interval "' output_ISI_name '"']}
    {'prog_path', '', ['--save-spike "' output_RAS_name '"']}
    {'save_conductance'}
    {'save_poisson_events'}
    {'extra_cmd', '%s', ''}
    % Recognized parameters that will not pass to the executable.
    {'neuron_model', '', ''}
    {'simu_method', '', ''}
    {'net', '', ''}
    {'net_adj', '', ''}
    {'cmd_str', '', ''}
    {'ext_T', '', ''}
    {'scee_mV', '', ''}
    {'scie_mV', '', ''}
    {'scei_mV', '', ''}
    {'scii_mV', '', ''}
    {'ps_mV', '', ''}
    {'psi_mV', '', ''}
    }.';
end

% Identify non-recognized parameters.
% Here slow. About 1.9ms in Octave, 0.8ms in Matlab 2014a.
opt_names = cellfun(@(x) x{1}, c_options, 'UniformOutput', false);
field_names = fieldnames(pm).';
for id_reg = find(ismember(field_names, opt_names) == false)
    warning('gen_neu:pm', 'pm.%s not recognized.', field_names{id_reg});
end

% Construct command line string.
c_options_str = cell(size(c_options));
cnt_options = 1;
for id_c = 1:length(c_options)
    pm_option = c_options{id_c};
    c_options_str{cnt_options} = get_option_str(pm, pm_option{:});
    if ~isempty(c_options_str{cnt_options})
        cnt_options = cnt_options + 1;
    end
end

pm.cmd_str = strjoin(c_options_str(1:cnt_options-1));
extra_data.cmdst = pm.cmd_str;
if mode_show_cmd
    disp(pm.cmd_str);
    return
end

% test if there is cached files.
have_data = exist(output_name, 'file') ...
         && exist(output_ISI_name, 'file') ...
         && exist(output_RAS_name, 'file');
if mode_extra_data
    have_data = have_data && exist(output_G_name, 'file') ...
         && exist(output_gating_name, 'file');
end
if (~have_data || new_run)...
   && ~mode_read_only...
   && (~mode_rm_only || mode_rm_only && nargout>0)
    % avoid data inconsistancy
    QuietDelete(output_ISI_name);
    QuietDelete(output_RAS_name);
    QuietDelete(output_G_name);
    QuietDelete(output_gating_name);
    if b_verbose
        fprintf('Parameter Preparing  : %.3f s\n', toc(t_start));
        fprintf('Command execution\n');
    end
    fflush(stdout);
    t_start = tic();
    rt = system(pm.cmd_str);           % Run simulation !
    if mode_run_in_background
        ISI=rt;
        return
    end
    if b_verbose
        fprintf('  Total              : %.3f s\n', toc(t_start));
    end
else
    if b_verbose
        fprintf('Parameter Preparing  : %.3f s\n', toc(t_start));
    end
    rt = 0;
end

t_start = tic;

if mode_read_only
    % Test if the file is exist and filled (modified more than 1 sec ago)?
    % Fixme: is there a way to test whether the files are using by
    %        other program in Matlab? lsof is not portable to M$ Windows
    % Note: dir().datenum only precise to seconds
    %       stat().mtime can precise to nano seconds, but not in Matlab.
    %f_info = stat(output_RAS_name);
    %if ~isempty(f_info) && (time() - f_info.mtime > 1.0)
    f_info = dir(output_RAS_name);
    if ~isempty(f_info) && (datenum(clock()) - f_info.datenum > 1.0/86400)
        rt = 0;      % we have data
    else
        return
    end
end

if rt ~= 0
    disp('pm.cmd_str =');
    disp(pm.cmd_str);
    error('Fail to generate data!');
end
% If requested, read and return data.
% At this point, the files should be generated, and ready for read.
if nargout > 0
    % Fixme: The Octave function isargout() is not in Matlab.
    %        So seems no way to skip non-required data in Matlab
    %        in some condition.
    if return_X_name
        X = output_name;
    else
        fid = fopen(output_name, 'r');
        X = fread(fid, [p, Inf], 'double');
        fclose(fid);
        len = size(X,2);
        if len ~= floor((pm.t+ext_T)/pm.stv)
            fprintf('gen_neu:\n');
            fprintf('  size(X,2) = %d (read in file),\n', size(X,2));
            fprintf('  floor((pm.t+ext_T)/pm.stv) = %d (expected)\n',...
                    floor((pm.t+ext_T)/pm.stv));
            warning('gen_neu:out', 'inconsistant data length! Exceed part will be truncated.');
        end
        len_cut = len - floor(pm.t/pm.stv);
        if (len_cut>0)
            X(:, 1:end-floor(pm.t/pm.stv)) = [];
        end
    end
    if (nargout>1)
        ISI = load('-ascii', output_ISI_name);
    end
    if (nargout>2)
        % We need to check if the file is non-empty (i.e. no any spike)
        % Otherwise load() may fail due to empty file.
        % TODO: use dir() in this part
        tmp_fd = fopen(output_RAS_name);
        rt1 = fseek(tmp_fd, 0, 'eof');
        if rt1 ~= 0
            fclose(tmp_fd);
            error('fail to call fseek(), check the file "%s"',...
                  output_RAS_name);
        end
        pos = ftell(tmp_fd);
        fclose(tmp_fd);
        if pos > 1
            ras = load('-ascii', output_RAS_name);
        else
            ras = [];
        end
        if exist('len_cut','var') && len_cut>0 && ~isempty(ras)
            ras(ras(:,2) <= len_cut*pm.stv, :) = [];
            ras(:,2) = ras(:,2) - len_cut*pm.stv;
        end
    end
    if mode_extra_data
        % read conductance and gating varialbes (if any)
        fid = fopen(output_G_name, 'r');
        extra_data.G = fread(fid, [2*p, Inf], 'double');
        extra_data.G(:, 1:end-floor(pm.t/pm.stv)) = [];
        fclose(fid);
        fid = fopen(output_gating_name, 'r');
        if fid >=0
          extra_data.gatings = fread(fid, [3*p, Inf], 'double');
          extra_data.gatings(:, 1:end-floor(pm.t/pm.stv)) = [];
          fclose(fid);
          % The order of gating variables is defined in
          %   single_neuron_dynamics.h
          % See id_h, id_m, id_n in struct Ty_HH_GH
        end
    end
end

% delete data if asked
if mode_rm_only
    QuietDelete(output_name);
    QuietDelete(output_ISI_name);
    QuietDelete(output_RAS_name);
    QuietDelete(output_G_name);
    QuietDelete(output_gating_name);
    if b_clean_net_file  % so pm.net was adjacency matrix
        % Remove the saved network file.
        QuietDelete(pm.net_path);
    end
end
QuietDelete(poisson_path);
QuietDelete(force_spike_path);

if b_verbose
    fprintf('Load data and clean  : %.3f s\n', toc(t_start));
end

end  % end of function gen_neu

% Give array a name.
% item = 'pr' 'ps' 'pri' or 'psi'.
function name = get_prps_str(pm, item)
    if has_nonempty_field(pm, item)
        if numel(pm.(item)) == 1
            name = [item '=' sprintf('%.5g', pm.(item))];
        else
            name = [item '=' BKDRHash(pm.(item))];
        end
    else
        name = '';
    end
end

function b = has_nonempty_field(stru, field_name)
    b = isfield(stru, field_name) && ~isempty(stru.(field_name));
end

% Convert struct field to arguement name and value in string.
% If only "field_name" is provided, format for contained data and option name
%   will auto determined.
%   Recognized formats include interger or floating scalar or vector, string.
%   If the field is an empty matrix,the option is treated as if not provided.
% If "formatting" or "option_name" is [], then it is treated as if not
%   provided. To pass empty string, use ''.
function str = get_option_str(s, field_name, formatting, option_name)
    if ~isfield(s, field_name) ||...
        isempty(s.(field_name)) && ~ischar(s.(field_name))
        str = '';
        return
    end
    % Auto determine the data formatting.
    if ~exist('formatting', 'var') || isempty(formatting) && ~ischar(formatting)
        f = s.(field_name);
        if isnumeric(f) && isvector(f)
            formatting = ' %.17g';
        elseif ischar(f)
            formatting = '"%s"';
        elseif isempty(f)
            formatting = '';
        else
            disp(f);
            error('Unhandled data type.');
        end
    end
    % Auto determine the command line option name.
    if ~exist('option_name', 'var') || isempty(option_name) && ~ischar(option_name)
        option_name = ['--' strrep(field_name , '_', '-')];
    end
    if isempty(formatting)
        % This 'if' is for damn matlab:
        % Error using sprintf
        % Function is not defined for sparse inputs.
        str1 = '';
    else
        str1 = sprintf(formatting, s.(field_name));
    end
    white_space = char(' '*ones(1, 1-(length(option_name)==0 || length(str1)==0 || str1(1)==' ')));
    str = [option_name white_space str1];
end

function z = inlineif(c, a, b)
    if c
        z = a;
    else
        if exist('b', 'var')
            z = b;
        else
            z = [];
        end
    end
end

%test
%{
  clear('pm');
  pm.net  = 'net_2_2';
  pm.neuron_model = 'LIF-GH';
  pm.simu_method = 'auto';
  pm.scee = 0.05;
  pm.ps   = 0.04;
  pm.pr   = 1.6;
  pm.t    = 1e4;
  pm.dt   = 1/32;
  pm.stv  = 0.5;
  [X, ISI, ras] = gen_neu(pm);
  % no error

  pm.nE = 1;
  [X, ISI, ras] = gen_neu(pm);
  % Number of neurons inconsist with the network!

  pm.nI = 2;
  [X, ISI, ras] = gen_neu(pm);

  pm = rmfield(pm, 'nI');
  pm = rmfield(pm, 'nE');
  X = gen_neu(pm, '');
  X = gen_neu(pm, '');
  X = gen_neu(pm, 'new');
  X = gen_neu(pm, 'rm');  size(X)
  X = gen_neu(pm, 'new,rm');  size(X)
  gen_neu(pm, '');
  gen_neu(pm, 'rm');

  gen_neu(pm,'');
  gen_neu(pm,'read');
  gen_neu(pm,'read,rm');

  X=gen_neu(pm,'nameX');
  X=gen_neu(pm,'nameX,rm');

  X=gen_neu(pm,'cmd');
  
  pm.net = [0 1; 0 0];
  pm.nE = 1;
  pm.nI = 1;
  pm.pr = [1.6 1.5];
  pm.extra_cmd = '--output-poisson-path poi.txt';
  X=gen_neu(pm, 'new,rm');
  poi = load('poi.txt');
  1/mean(diff(poi(poi(:,1)==1, 2)))
  1/mean(diff(poi(poi(:,1)==2, 2)))
  
  pm.ps = [0.04 0.05];
  X=gen_neu(pm, 'new,rm');  % no error
  poi = load('poi.txt');
  1/mean(diff(poi(poi(:,1)==1, 2)))
  1/mean(diff(poi(poi(:,1)==2, 2)))
  all(poi(poi(:,1)==1, 3) == pm.ps(1))
  all(poi(poi(:,1)==2, 3) == pm.ps(2))
%}

% vim: et sw=4 sts=4
