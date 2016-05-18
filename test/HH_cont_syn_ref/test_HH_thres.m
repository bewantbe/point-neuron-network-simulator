% Single neuron case
% test HH model spike

run('/home/xyy/matcode/GC_clean/startup.m');
addpath('/home/xyy/matcode/prj_GC_clean/HH/');

addpath('/home/xyy/code/point-neuron-network-simulator/mfile/');

    is_octave = exist('OCTAVE_VERSION','builtin') ~= 0;
    if is_octave
        fontsize = 22;
        linewidth = 3;
    else
        fontsize = 20;
        linewidth = 1;
    end
    pic_prefix = 'pic_tmp/';
    pic_output_bw = @(st)print('-deps',  [pic_prefix, st, '.eps']);
    pic_output    = @(st)print('-depsc2',[pic_prefix, st, '.eps']);
%    set(0, 'defaultfigurevisible', 'off');
    set(0, 'defaultlinelinewidth', linewidth);
    set(0, 'defaultaxesfontsize', fontsize);

% Get reference answer
clear('pm');
pm.neuron_model = '/home/xyy/code/point-neuron-network-simulator/test/HH_cont_syn_ref/IFsimu_HH/bin/Release/raster_tuning_HH3_gcc';
%pm.net  = 'net_1_0';
pm.net  = 'net_1_0';
pm.nI   = 0;
pm.scee = 0.05;
pm.pr   = 1.0;
pm.ps   = 0.02;
pm.t    = 1e4;
pm.dt   = 1.0/32;
pm.stv  = pm.dt;
pm.seed = 21414;
pm.extra_cmd = '-q';
thres = 84.0;
seting_st = sprintf('HH3_pr=%.1g_ps=%.1g_thres=%.0f', pm.pr, pm.ps, thres);
tic
[X_ref, ISI_ref, ras_ref] = gen_HH(pm, 'rm');
toc
X_ref(:, 1) = [];  % somehow the first data point is incorrect.
X_ref = 10*(X_ref-6.5);
1e3 ./ ISI_ref

id_neu = 1;

id_from = id_neu;
id_to = id_neu;
t_range = [-1.5,5];
[Y, T, id_s] = SpikeSlicing(X_ref, ras_ref, pm.stv, id_to, id_from, t_range);
figure(45);
plot(T', Y');
xlabel('t_{rel} /ms');
ylabel('volt (mV)');
pic_output(sprintf('spikes_%s', seting_st));

Z = X_ref(id_to, :);
id_subthres = true(size(Z));
id_subthres(id_s') = false;
Z = Z(id_subthres);
Z(:, end-10/pm.stv:end) = [];
Z(:, 1:10/pm.stv) = [];
figure(30);
plot(Z);
xlabel('subthreshold time point');
ylabel('volt (mV)');
pic_output(sprintf('sub-volt_%s', seting_st));

figure(31);
hist(Z, 100);
ylim([0 100]);
xlabel('volt (mV)');
ylabel('count');
pic_output(sprintf('volt-hist_%s', seting_st));

