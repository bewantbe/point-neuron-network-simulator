% Pick spike events in neuron set id_neu and time range t_range.
% id_neu list the neuron index wantted in a row.
%        if id_neu has second row, further rewrite neuron index to it.
% t_range: length=1: select range [0 t_range]
%          length=2: select range [t_range(1) t_range(2)] and subtract
%                    t_range(1) from the timings.
%          length=3: select range [t_range(1) t_range(2)] and subtract
%                    t_range(3) from the timings.
% Ex:
%    ras = [randi(5, 20, 1), 50*rand(20, 1)];
%    ras(:, [2 1]) = sortrows(ras(:, [2 1]));   % has to be time sorted
%    ras312 = ras_pick(ras, [1 2 3; 3 1 2], ras(end, 2)/2);

function local_ras = ras_pick(ras, id_neu, t_range)
  if isempty(ras)
    return
  end

  % Default t_range.
  if ~exist('t_range', 'var') || isempty(t_range)
    t_range = ras(end, 2);
  end
  if numel(t_range) == 1
    t_range = [0 t_range];
  end
  if numel(t_range) == 2
    t_range(3) = t_range(1);
  end

  % Select time range.
  if ras(1, 2) < t_range(1) || t_range(2) < ras(end, 2)
    % speed: 0.07 sec / 1e7
    id_bg = find(ras(:,2) >= t_range(1), 1, 'first');
    id_ed = find(ras(:,2) <= t_range(2), 1, 'last');
    local_ras = ras(id_bg:id_ed, :);
  else
    local_ras = ras;
  end

  p = max(ras(:,1));  % speed: 0.02 sec / 1e7

  % Default id_neu.
  if ~exist('id_neu') || isempty(id_neu)
    id_neu = 1:p;
  elseif any(id_neu(1, :) > p)
    % No such index, probably due to silent neuron.
    % Without this id_b_neu and id_map will be auto expanded.
    id_neu(:, id_neu(1, :) > p) = [];
  end
  
  % Select variables.
  id_b_neu = false(1, p);
  id_b_neu(id_neu(1, :)) = true;
  if ~all(id_b_neu)
    % speed: 0.05~0.18 sec / 1e7
    local_ras = local_ras(id_b_neu(local_ras(:,1)), :);
  end
  
  % Rewrite ras index.
  if size(id_neu, 1) >= 2
    % speed: 0.12 sec / 1e7
    id_map = zeros(p, 1);
    id_map(id_neu(1, :)) = id_neu(2, :);
    local_ras(:, 1) = id_map(local_ras(:, 1));
  end
  
  % Rewrite time.
  if t_range(3) ~= 0
    local_ras(:, 2) = local_ras(:, 2) - t_range(3);
  end
end
