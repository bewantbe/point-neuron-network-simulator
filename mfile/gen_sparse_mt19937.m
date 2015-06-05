
function net = gen_sparse_mt19937(p, sparseness, seed)
  s = randMT19937('state');
  if isstruct(p)
    sparseness = p.sparseness;
    seed       = p.seed;
    p          = p.p;
  end
  randMT19937('state', seed);
  net = randMT19937(p) < sparseness;
  net(eye(p)==1) = 0;
  randMT19937('state', s);
