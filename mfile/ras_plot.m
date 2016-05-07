function ras_plot(ras, t_bg, t_ed, id_neu)
if ~exist('t_bg', 'var') || isempty(t_bg)
  t_bg = 0;
end
if ~exist('t_ed', 'var') || isempty(t_ed)
  t_ed = ras(end, 2);
end
p = max(ras(:,1));
if ~exist('id_neu', 'var') || isempty(id_neu)
  id_neu = 1:p;
end

% select time range
id_bg = find(ras(:,2) >= t_bg, 1, 'first');
id_ed = find(ras(:,2) <= t_ed, 1, 'last');
local_ras = ras(id_bg:id_ed, :);

% select variable
id_b_neu = false(1, p);
id_b_neu(id_neu) = true;
local_ras = local_ras(id_b_neu(local_ras(:,1)), :);

t_len = t_ed - t_bg;
cla
line([local_ras(:, 2)'; local_ras(:, 2)']-t_bg,...
     [local_ras(:, 1)'; local_ras(:, 1)'-0.7], 'color', [0 0 0], 'linewidth', 5);
xlabel('time (ms)');
ylabel('neuron id');

axis([0, t_len, min(id_neu)-1, max(id_neu)]);

