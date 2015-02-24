
# Generate test (answer) data
make -f makefile_gcc
% case 1
/home/xyy/code/IFsimu_hg/raster_tuning -inf - -ng -t 1e5 -dt 0.5 -stv 0.5 -n 2 -mat - -scee 0.01 -scei 0.01 -scie 0.02 -scii 0.02 -pr 1 -ps 0.012 --RC-filter 0 1 --bin-save -o a_volt.dat --save-spike a_ras.txt
% case 2
/home/xyy/code/IFsimu_hg/raster_tuning -inf - -ng -t 1e4 -dt 0.5 -stv 0.5 -n 2 2 -mat - -scee 0.01 -scei 0.012 -scie 0.018 -scii 0.02 -pr 1 -ps 0.011 --RC-filter 0 1 --bin-save -o a_volt.dat --save-spike a_ras.txt
% case 3
/home/xyy/code/IFsimu_hg/raster_tuning -inf - -ng -t 1e3 -dt 0.5 -stv 0.5 -n 100 -mat net_100_0X65763652.txt -scee 0.005 -scei 0.005 -scie 0.01 -scii 0.01 -pr 1 -ps 0.009 --RC-filter 0 1 --bin-save -o a_volt.dat --save-spike a_ras.txt
% case 4
time /media/tsb/home/xyya/matcode/prj_GC_clean/HH/raster_tuning_LIF -inf - -ng -t 1e3 -dt 0.5 -stv 0.5 -n 700 300 -mat - -scee 0.002 -scei 0.002 -scie 0.004 -scii 0.004 -pr 1 -ps 0.005 --RC-filter 0 1 --bin-save -o a_volt.dat --save-spike a_ras.txt
%real	0m32.966s

# Generate data to be verified
time g++ -g -O0 -std=c++11 -Wall main.cpp math_helper.cpp -lboost_program_options -o bin/vec_IFsimu
% case 1
bin/vec_IFsimu --t 1e5 --dt 0.5 --stv 0.5 --nE 2 --net - --scee 0.01 --scie 0.01 --scei 0.02 --scii 0.02 --ps 0.012 --initial-state-path /home/xyy/code/vec_IFsimu/neu_state_init.txt --input-event-path /home/xyy/code/vec_IFsimu/poisson_events.txt -o v_volt.dat --conductance-path=v_cond.dat --ras-path=v_ras.txt
% case 2
bin/vec_IFsimu --t 1e4 --dt 0.5 --stv 0.5 --nE 2 --nI 2 --net - --scee 0.01 --scie 0.012 --scei 0.018 --scii 0.02 --ps 0.011 --initial-state-path /home/xyy/code/vec_IFsimu/neu_state_init.txt --input-event-path /home/xyy/code/vec_IFsimu/poisson_events.txt -o v_volt.dat --conductance-path=v_cond.dat --ras-path=v_ras.txt
% case 3
bin/vec_IFsimu --t 1e3 --dt 0.5 --stv 0.5 --nE 100 --net net_100_0X65763652.txt --scee 0.005 --scie 0.005 --scei 0.01 --scii 0.01 --ps 0.009 --initial-state-path /home/xyy/code/vec_IFsimu/neu_state_init.txt --input-event-path /home/xyy/code/vec_IFsimu/poisson_events.txt -o v_volt.dat --conductance-path=v_cond.dat --ras-path=v_ras.txt
% case 4
bin/vec_IFsimu --t 1e3 --dt 0.5 --stv 0.5 --nE 700 --nI 300 --net - --scee 0.002 --scie 0.002 --scei 0.004 --scii 0.004 --pr 1 --ps 0.005 -o v_volt.dat --ras-path=v_ras.txt --isi-path=v_isi.txt
%real	0m42.141s
%real	0m20.069s
%real	0m18.427s
%real	0m15.403s
%real	0m14.986s  % hand tune
%real	0m14.516s  % Lagrange form
%real	0m14.339s  % poly coeff
%real	0m9.692s   % 6 binary + 3 newton's

# Compare
X0 = load('neu_state_init.txt');
p = size(X0,2)/3;
fid = fopen('v_volt.dat');
X1 = fread(fid, [p Inf], 'double');
fclose(fid);
fid = fopen('v_cond.dat');
c_g = fread(fid, [p*2 Inf], 'double');
fclose(fid);


a_v = X0(1:end, [3*(0:p-1)+1])';
a_gE = X0(1:end, [3*(0:p-1)+2])';
a_gI = X0(1:end, [3*(0:p-1)+3])';
c_v = X1;
c_gE = c_g(1:2:end, :);
c_gI = c_g(2:2:end, :);

rg = (1:190);
rg_p = 1:p;
figure(1);
plot(rg, a_v(rg_p,rg), rg, c_v(rg_p,rg));

figure(1);
plot(rg, a_v(rg_p,rg) - c_v(rg_p,rg));


figure(5);
plot(rg, a_gE(:,rg), '-o');

figure(6);
plot(rg, a_gE(:,rg), '-o');

figure(3);
plot(rg, a_gE(rg_p,rg), rg, c_gE(rg_p,rg));


rg = (1:200);
rg_p = 1:p;
figure(3);
plot(rg, a_gE(rg_p,rg) - c_gE(rg_p,rg));
%legend('1','2','3','4');

% RAS
a_ras = load('a_ras.txt');
v_ras = load('v_ras.txt');
figure(8);
plot(a_ras(1:5000,1) == v_ras(1:5000,1));

figure(10);
clf;
cla;
%local_ras = ras(1:1000, :);
local_ras = ras;
line([local_ras(:, 2)'; local_ras(:, 2)'],...
     [local_ras(:, 1)'; local_ras(:, 1)'-0.8]);

figure(11);
clf;
cla;
local_ras = ras(1:1000, :);
line([local_ras(:, 2)'; local_ras(:, 2)'],...
     [local_ras(:, 1)'; local_ras(:, 1)'-0.8]);

