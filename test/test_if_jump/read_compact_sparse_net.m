% Read compact sparse format

function w = read_compact_sparse_net(net_path)

nv=load(net_path);
p = find(nv - (length(nv)-1:-1:0)' == 0, 1) - 1;
cntj = nv(1:p+1);
wj = nv(p+2:end) + 1;
wi = zeros(size(wj));
for k = 1 : length(cntj)-1
  wi(cntj(k)+1:cntj(k+1)) = k;
end
w = sparse(wj, wi, 1);
w(p,p) = 0;

