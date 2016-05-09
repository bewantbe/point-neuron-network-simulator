%

%[X_ref, ISI_ref, ras_ref] = gen_neu(pm, 'new,rm');
%t_range = 20;
%id_from = 1;
%id_to = 1;
%[Y, T] = SpikeSlicing(X_ref, ras_ref, pm.stv, id_to, id_from, t_range);
%plot(T', Y');

function [Y, T] = SpikeSlicing(X, ras, stv, id_to, id_from, t_range)

t_end = size(X,2)*stv;
t_c = ras(ras(:,1)==id_from, 2);  % time of id_from spikes

if numel(t_range)==1
  t_l = -t_range/2;
  t_u =  t_range/2;
else  % numel(t_range)==2
  t_l = t_range(1);
  t_u = t_range(2);
end

% remove boundary
t_c( t_c + t_l <= 0.5*stv | t_c + t_u >= t_end - 0.5*stv ) = [];

n_s = floor((t_u - t_l)/stv);
n_l = fix(t_l/stv);

i_c = round(t_c / stv);

id_s = bsxfun(@plus, n_l + (0:n_s-1), i_c);

Y = reshape(X(id_to, id_s), size(id_s));
T = bsxfun(@plus, id_s*stv, -t_c);

