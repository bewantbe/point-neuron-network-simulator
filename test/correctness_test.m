%

cmd_correctness_prog = '/home/xyy/code/IFsimu_hg/raster_tuning';
cmd_speed_ref_prog = '/media/tsb/home/xyya/matcode/prj_GC_clean/HH/raster_tuning_LIF';

s_cmd_correctness_ref = {
'-inf - -ng -t 1e5 -dt 0.5 -stv 0.5 -n 2 -mat - -scee 0.01 -scei 0.01 -scie 0.02 -scii 0.02 -pr 1 -ps 0.012 --RC-filter 0 1 --bin-save -o a_volt.dat --save-spike a_ras.txt'
'-inf - -ng -t 1e4 -dt 0.5 -stv 0.5 -n 2 2 -mat - -scee 0.01 -scei 0.012 -scie 0.018 -scii 0.02 -pr 1 -ps 0.011 --RC-filter 0 1 --bin-save -o a_volt.dat --save-spike a_ras.txt'
'-inf - -ng -t 1e3 -dt 0.5 -stv 0.5 -n 100 -mat net_100_0X65763652.txt -scee 0.005 -scei 0.005 -scie 0.01 -scii 0.01 -pr 1 -ps 0.009 --RC-filter 0 1 --bin-save -o a_volt.dat --save-spike a_ras.txt'
};
s_cmd_speed_ref = {
'-inf - -ng -t 1e3 -dt 0.5 -stv 0.5 -n 700 300 -mat - -scee 0.002 -scei 0.002 -scie 0.004 -scii 0.004 -pr 1 -ps 0.005 --RC-filter 0 1 --bin-save -o a_volt.dat --save-spike a_ras.txt'
};

cmd_correctness_target_prog = '../bin/gen_neu';
cmd_speed_target_prog = '../bin/gen_neu';

s_cmd_correctness_target  = {
'--t 1e5 --dt 0.5 --stv 0.5 --nE 2 --net - --scee 0.01 --scie 0.01 --scei 0.02 --scii 0.02 --ps 0.012 --initial-state-path /home/xyy/code/vec_IFsimu/neu_state_init.txt --input-event-path /home/xyy/code/vec_IFsimu/poisson_events.txt -o v_volt.dat --conductance-path=v_cond.dat --ras-path=v_ras.txt'
'--t 1e4 --dt 0.5 --stv 0.5 --nE 2 --nI 2 --net - --scee 0.01 --scie 0.012 --scei 0.018 --scii 0.02 --ps 0.011 --initial-state-path /home/xyy/code/vec_IFsimu/neu_state_init.txt --input-event-path /home/xyy/code/vec_IFsimu/poisson_events.txt -o v_volt.dat --conductance-path=v_cond.dat --ras-path=v_ras.txt'
'--t 1e3 --dt 0.5 --stv 0.5 --nE 100 --net net_100_0X65763652.txt --scee 0.005 --scie 0.005 --scei 0.01 --scii 0.01 --ps 0.009 --initial-state-path /home/xyy/code/vec_IFsimu/neu_state_init.txt --input-event-path /home/xyy/code/vec_IFsimu/poisson_events.txt -o v_volt.dat --conductance-path=v_cond.dat --ras-path=v_ras.txt'
};

s_cmd_speed_target = {
'--t 1e3 --dt 0.5 --stv 0.5 --nE 700 --nI 300 --net - --scee 0.002 --scie 0.002 --scei 0.004 --scii 0.004 --pr 1 --ps 0.005 -o v_volt.dat --ras-path=v_ras.txt --isi-path=v_isi.txt'
'--neuron-model LIF-GH --t 1e3 --dt 0.125 --stv 0.5 --nE 700 --nI 300 --net - --scee 0.004 --scie 0.004 --scei 0.008 --scii 0.008 --pr 1 --ps 0.020 -o v_volt.dat --ras-path=v_ras.txt --isi-path=v_isi.txt'
'--neuron-model HH-GH --t 1e3 --dt 0.03125 --stv 0.5 --nE 700 --nI 300 --net - --scee 0.006 --scie 0.006 --scei 0.012 --scii 0.012 --pr 1 --ps 0.030 -o v_volt.dat --ras-path=v_ras.txt --isi-path=v_isi.txt'
};

%'--neuron-model HH-GH --t 1.48e3 --dt 0.03125 --stv 0.5 --nE 1 --nI 0 --net - --scee 0.000 --scie 0.000 --scei 0.000 --scii 0.000 --pr 1 --ps 0.040 -o v_volt.dat --ras-path=v_ras.txt --isi-path=v_isi.txt'

%system(['time ' cmd_speed_target_prog ' --neuron-model LIF-G ' s_cmd_speed_target{1}]);

%system(['time ' cmd_speed_target_prog ' --neuron-model LIF-G-Sparse ' s_cmd_speed_target{1}]);

pamcomm = '--t 1e4 --dt 0.5 --stv 0.5 --nE 200 --net - --scee 0.0 --scie 0.0 --scei 0.0 --scii 0.0 --ps 0.012';

%pamcomm = '--t 1e3 --dt 0.5 --stv 0.5 --nE 2 --net - --scee 0.01 --scie 0.01 --scei 0.02 --scii 0.02 --ps 0.012';

pam0 = ' -o r_volt.dat --conductance-path=r_cond.dat --ras-path=r_ras.txt';
system(['time ' cmd_speed_target_prog ' --neuron-model LIF-G ' pamcomm pam0]);
pam1 = ' -o v_volt.dat --conductance-path=v_cond.dat --ras-path=v_ras.txt';
system(['time ' cmd_speed_target_prog ' --neuron-model LIF-G-Sparse ' pamcomm pam1]);

p = 200;
stv = 0.5;
len = 1e4/stv;

s_t = stv*(1:len);
X0 = ReadDouble('r_volt.dat', p);
X1 = ReadDouble('v_volt.dat', p);
ras0 = load('r_ras.txt');
ras1 = load('v_ras.txt');

ST0 = SpikeTrains(ras0, p, len, stv);
ST1 = SpikeTrains(ras1, p, len, stv);

aX0 = X0 + cumsum(ST0, 2);
aX1 = X1 + cumsum(ST1, 2);

figure(1);
plot(s_t, aX0-aX1);
title('X0-X1');

