% Generator for Erdős–Rényi random graphs
% https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model

function net = gen_net_er(p, sparseness, seed)
  if isstruct(p)
    sparseness = p.sparseness;
    seed       = p.seed;
    p          = p.p;
  end
  s = randMT19937('state');
  randMT19937('state', seed);
  net = randMT19937(p) < sparseness;
  randMT19937('state', s);
  net(eye(p)==1) = 0;

