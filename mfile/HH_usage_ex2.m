% Usage example.
% Single neuron simulation: With sine current input.

pm = [];
pm.neuron_model = 'HH-GH-sine';
pm.simu_method  = 'simple';
pm.net   = 'net_1_0';
pm.pr    = 0;
pm.ps    = 0;
pm.t     = 1e3;
pm.dt    = 1.0/32;
pm.stv   = 0.125;
pm.seed  = 235478;
pm.extra_cmd = '';
pm.sine_amp  = 1;

s_amp  = linspace(1, 2, 4);
s_freq = linspace(0, 0.25, 40);
s_v_max = zeros(length(s_freq), length(s_amp));

tic
for id_amp = 1:length(s_amp)
  for id_freq = 1:length(s_freq)
    pm.sine_amp  = s_amp(id_amp);
    pm.sine_freq = s_freq(id_freq);
    X = gen_neu(pm, 'new,rm,ext_T');
    s_v_max(id_freq, id_amp) = max(X);
  end
end
toc

fontsize = 22;
linewidth = 3;
set(0, 'defaultlinelinewidth', linewidth);
set(0, 'defaultaxesfontsize', fontsize);

figure(1);
plot(s_freq*1000, s_v_max, '-o');
ylabel('V_{max}');
xlabel('Hz');
title('frequency resonance of V');

amp_big = 5;
%s_freq = linspace(0.127, 0.140, 1000);  % Chaos ?
c_ISI = cell(1, length(s_freq));

tic
for id_freq = 1:length(s_freq)
  pm.sine_amp  = amp_big;
  pm.sine_freq = s_freq(id_freq);
  [~, ~, ras] = gen_neu(pm, 'new,rm,ext_T');
  if id_amp == length(s_amp)
    if size(ras, 1) >= 2
      c_ISI{id_freq} = diff(ras(:, 2));
    else
      c_ISI{id_freq} = [];
    end
  end
end
toc

figure(3);
hold off
cla
clf
hold on
for k=1:length(c_ISI)
  if ~isempty(c_ISI{k})
    scatter(1000*s_freq(k)*ones(size(c_ISI{k})), c_ISI{k}, 10);
  end
end
for k = 1:6
  plot(1000 * s_freq, k ./ s_freq, 'color', [0.9 0.9 0.9]);
end
ylim([0 100]);
ylabel('ISI (ms)');
xlabel('Input Freq (Hz)');
title('frequency resonance of ISI');

