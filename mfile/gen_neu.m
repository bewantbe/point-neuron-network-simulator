% Generate HH neuron data by calling gen_neu
% Will use cached data automatically
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
%  'read'   Read data files if exist, otherwise do nothing and return [];
%  'nameX'  return path to the voltage data file, instead of read it;
%  'cmd'    Show command call to raster_tuning_HH, then exit. Useful for debug
%  'ext_T'  Generate a bit more data, so reduce "head effect".
%  'v'      Be verbose, with 'cmd' will also run and show cmd.
%  'h'      Show help for the command line generator.
%
% return value pm is a "normalized" input parameter set

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
    error('Lack of input parameters.');
end
if ~exist('fflush', 'builtin')
    fflush = @(a) 0;
    stdout = 0;
end
X=[];
ISI=[];
ras=[];

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
        if isfield(pm, [fv '_mV'])
            s = pm.([fv '_mV']) * PSP_v(id_fv);
            if isfield(pm, fv) && s ~= pm.(fv)
                error(['incompatible EPSP/IPSP strength specification: ' fv ' ' fv '_mV']);
            end
            pm.(fv) = s;
        end
    end
end

if ~ischar(pm.net)
  b_clean_net_file = true;
else
  b_clean_net_file = false;
end

if ~isfield(pm, 'prog_path')
    % Search gen_neu in "../bin" then "./"
    pathdir = fileparts(mfilename('fullpath'));
    exepath = sprintf('%s%s..%sbin%sgen_neu', pathdir, filesep, filesep, filesep);
    if ~exist(exepath, 'file')
        exepath = sprintf('%s%sgen_neu', pathdir, filesep);
    end
else
    exepath = pm.prog_path;
end

% Default generator settings
new_run        = false;
return_X_name  = false;
mode_rm_only   = false;
mode_show_cmd  = false;
mode_read_only = false;
mode_extra_data = false;
mode_run_in_background = false;

if ~isfield(pm, 'extra_cmd')
    pm.extra_cmd = '';
end
pm.extra_cmd = [' ' strtrim(pm.extra_cmd) ' '];
b_verbose = ~isempty([strfind(pm.extra_cmd,' -v ') strfind(pm.extra_cmd,' --verbose ')]);  % show time cost
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
        system(['"' exepath '" -h']);
        return
    otherwise
        error('no this option: "%s"', tok);
    end
end
pm.ext_T = ext_T;

% Default parameter values
if ~exist('data_dir_prefix', 'var')
    data_dir_prefix = ['.', filesep, 'data', filesep];
end
if ~isfield(pm, 'net') || isempty(pm.net)
    pm.net = 'net_1_0';
end
if ischar(pm.net)
    [network, mat_path] = get_network(pm.net, data_dir_prefix);
    if isempty(network)
        if isfield(pm, 'net_adj') && isnumeric(pm.net_adj) ...
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
if ~isfield(pm, 'nI') || isempty(pm.nI)
    pm.nI = 0;  % number of inhibitory neurons
end
if isfield(pm, 'nE') && ~isempty(pm.nE)
    if pm.nI + pm.nE ~= p
        fprintf('  pm.nI + pm.nE = %d, p_net = %d\n', pm.nI + pm.nE, p);
        error('gen_neu: Number of neurons inconsist with the network!');
    end
else
    pm.nE = p - pm.nI;
end
if ~isfield(pm, 'stv') || isempty(pm.stv)
    pm.stv = 0.5;
end
if ~isfield(pm, 'dt') || isempty(pm.dt)
    pm.dt = 1.0/32;
end
if ~isfield(pm, 'scee') || isempty(pm.scee)
    pm.scee = 0;
end
if ~isfield(pm, 'scie') || isempty(pm.scie)
    pm.scie = 0;  % Strength from Ex. to In.
end
if ~isfield(pm, 'scei') || isempty(pm.scei)
    pm.scei = 0;  % Strength from In. to Ex.
end
if ~isfield(pm, 'scii') || isempty(pm.scii)
    pm.scii = 0;
end
s_tmp = strtrim(pm.extra_cmd);
if ~isempty(s_tmp) && strcmp(s_tmp(end), '&') == 1
    % start the data generation in background, then return immediately
    mode_run_in_background = true;  
end
if ~isfield(pm, 'neuron_model') || isempty(pm.neuron_model)
    error('neuron_model not specified!');
end
if ~isfield(pm, 'simu_method') || isempty(pm.simu_method)
    disp('Warning: .simu_method not set! Using auto mode.');
    pm.simu_method = 'auto';
end
neuron_model_name = pm.neuron_model;
program_name = sprintf(...
    '"%s" --neuron-model %s --simulation-method %s',...
    exepath, pm.neuron_model, pm.simu_method);

if ~isfield(pm, 'st_extra_inf_post')
    pm.st_extra_inf_post = '';
end

if xor(isfield(pm,'pr'), isfield(pm,'ps'))
    warning('pr and ps not privided together.');
end
if xor(isfield(pm,'pri'), isfield(pm,'psi'))
    warning('pri and psi not privided together.');
end

[pr_str pr_name] = get_prps_str(pm, 'pr');
[ps_str ps_name] = get_prps_str(pm, 'ps');
[pri_str pri_name] = get_prps_str(pm, 'pri');
[psi_str psi_name] = get_prps_str(pm, 'psi');

% construct file paths
st_sc = strrep(mat2str([pm.scee, pm.scie, pm.scei, pm.scii]),' ',',');
st_p  = strrep(mat2str([pm.nE, pm.nI]),' ',',');
file_inf_st =...
    sprintf('%s_p=%s_sc=%s_%s_%s_%s_%s_stv=%g_t=%.2e',...
            pm.net, st_p(2:end-1), st_sc(2:end-1), pr_name,...
            ps_name, pri_name, psi_name, pm.stv, pm.t + ext_T);
file_prefix = [data_dir_prefix, neuron_model_name, '_'];
output_name     = [file_prefix, 'volt_',file_inf_st,'.dat'];
output_ISI_name = [file_prefix, 'ISI_', file_inf_st,'.txt'];
output_RAS_name = [file_prefix, 'RAS_', file_inf_st,'.txt'];
output_G_name   = [file_prefix, 'G_',file_inf_st,'.dat'];
output_gating_name = [file_prefix, 'gating_',file_inf_st,'.dat'];

% construct command string
% ! NOTE: pm.scie is strength from Ex. to In. type
%         Seems biologist love this convention, I have no idea.
st_neu_s =...
    sprintf('--scee %.16e --scie %.16e --scei %.16e --scii %.16e',...
            pm.scee, pm.scie, pm.scei, pm.scii);
st_neu_param =...
    sprintf('--nE %d --nI %d --net "%s" %s %s %s %s %s',...
            pm.nE, pm.nI, mat_path, pr_str, ps_str, pri_str, psi_str, st_neu_s);
if issparse(pm.net_adj)
  st_neu_param = [st_neu_param ' --sparse-net'];
end
if isfield(pm, 'sine_amp')
    st_neu_param = [st_neu_param,...
      sprintf(' --current-sine-amp %.16e', pm.sine_amp)];
end
if isfield(pm, 'sine_freq')
    st_neu_param = [st_neu_param,...
      sprintf(' --current-sine-freq %.16e', pm.sine_freq)];
end
if isfield(pm, 'synaptic_delay') && ~isempty(pm.synaptic_delay)
    st_neu_param = [st_neu_param,...
      sprintf(' --synaptic-delay %.16e', pm.synaptic_delay)];
end
if isfield(pm, 'synaptic_net_delay') && ~isempty(pm.synaptic_net_delay)
    net_delay_path = save_network(pm.synaptic_net_delay, ...
                     [data_dir_prefix 'net_delay_']);
    st_neu_param = [st_neu_param,...
      sprintf(' --synaptic-net-delay %s', net_delay_path)];
end
if isfield(pm, 'extI')
    st_neu_param = [st_neu_param,...
      sprintf(' --extI %.16e', pm.extI)];
end

if isfield(pm, 'input_event')
    [~, tmp_f_name] = fileparts(tempname('./'));
    poisson_path = [data_dir_prefix 'external_event_' tmp_f_name '.txt'];
    pm.extra_cmd = [pm.extra_cmd ' --input-event-path ' poisson_path];
    if (size(pm.input_event,2) ~= 2 && size(pm.input_event,2) ~= 3)
        error('pm.input_event must be 2 or 3 column, the order is "id time [strength]"');
    end
    func_fraction = @(x) x - floor(x);
    if any(func_fraction(pm.input_event(:,1)))
        error('pm.input_event first column must be neuron index(0 based).');
    end
    if ~issorted(pm.input_event(:,2))
        warning('pm.input_event is not time sorted, sorting for you.');
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
    poisson_path = [];
end

if isfield(pm, 'force_spikes')
    [~, tmp_f_name] = fileparts(tempname(['.' filesep]));
    force_spike_path = [data_dir_prefix tmp_f_name 'force_spike.txt'];
    pm.extra_cmd = [pm.extra_cmd ' --force-spike-list ' force_spike_path];
    fd = fopen(force_spike_path, 'w');
    fprintf(fd, '%d %.6f\n', pm.force_spikes');
    fclose(fd);
else
    force_spike_path = [];
end

st_sim_param =...
    sprintf('--t %.16e --dt %.17e --stv %.17e',...
            pm.t + ext_T, pm.dt, pm.stv);
if isfield(pm, 'seed') && ~isempty(pm.seed) && strcmpi(pm.seed, 'auto')==0
    if isnumeric(pm.seed)
        str = sprintf(' %lu', pm.seed);
    else
        str = pm.seed;
    end
    st_sim_param = [st_sim_param, sprintf(' --seed %s', str)];
else
    st_sim_param = [st_sim_param, ' --seed-auto'];
end
st_paths =...
    sprintf(' --volt-path "%s" --isi-path "%s" --ras-path "%s"',...
            output_name, output_ISI_name, output_RAS_name);
if mode_extra_data
    st_paths = [st_paths ...
        sprintf(' --conductance-path "%s" --ion-gate-path "%s"',...
                output_G_name, output_gating_name)];
end
cmdst = sprintf('%s %s %s %s %s',...
                program_name, st_neu_param, st_sim_param, st_paths, pm.extra_cmd);
extra_data.cmdst = cmdst;
pm.cmd_str = cmdst;
if mode_show_cmd
    disp(cmdst);
    return
end

% test if there is cached files.
have_data = exist(output_name, 'file') ...
         && exist(output_ISI_name, 'file') ...
         && exist(output_RAS_name, 'file');
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
    rt = system(cmdst);           % Run simulation !
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
    error('Fail to generate data!');
end
% If required, read and return data
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
            warning('inconsistant data length! Exceed part will be truncated.');
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
if ~isempty(poisson_path)
    QuietDelete(poisson_path);
end
if ~isempty(force_spike_path)
    QuietDelete(force_spike_path);
end

if b_verbose
    fprintf('Load data and clean  : %.3f s\n', toc(t_start));
end

end  % end of function

% Convert array to command line option
% item = 'pr' 'ps' 'pri' or 'psi'.
function [str name] = get_prps_str(pm, item)
    if isfield(pm, item) && ~isempty(pm.(item))
        str = ['--' item];
        str = [str, sprintf(' %.17g', pm.(item))];
        if numel(pm.(item)) == 1
            name = [item '=' sprintf('%g', pm.(item))];
        else
            name = [item '=' BKDRHash(pm.(item))];
        end
    else
        str = '';
        name = '';
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
