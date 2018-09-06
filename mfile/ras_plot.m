% Plot raster (spiking events).
% Usage:
%   hd = ras_plot(ras, true);
%   set(hd, 'linewidth', 2);
%   xlabel('time (ms)');
%   ylabel('neuron id');
% To restore the old behaviour:
%   hd = ras_plot(ras, t_bg, t_ed, id_neu, y_scale);
% Write this:
%   hd = ras_plot( ras_pick(ras, id_neu, [t_bg, t_ed]) );

function hd = ras_plot(ras, ras_line_mode)
if isempty(ras)
  cla;
end
if ~exist('ras_line_mode', 'var') || isempty(ras_line_mode)
  ras_line_mode = max(ras(:,1)) <= 100;
end

if ras_line_mode
  cla;
  hd = ...
  line([ras(:, 2)'; ras(:, 2)'],...
       [ras(:, 1)'; ras(:, 1)'-0.7], ...
       'color', [0 0 0]);
else
  hd = plot(ras(:, 2), ras(:, 1), '.');
end

