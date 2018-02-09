% Save the adjacency matrix `A' in path `path_prefix' with a hashed file name:
%   matpath = save_network(A, path_prefix);
% To read the matrix from the file:
%   A = load('-ascii', matpath);

function [matpath, matname] = save_network(A, path_prefix)

% sane test of the matrix A
if ~exist('A','var') || isempty(A) || (~isnumeric(A) && ~islogical(A))
  error('Input a numeric matrix please!');
end
if ~ismatrix(A) || size(A,1) ~= size(A,2)
  error('Not a square matrix!');
end
A = double(A);  % int32,logical to double, so get consistent hash result

if ~exist('path_prefix','var')
  path_prefix = '';    % save to default dir (usually working dir)
end
pathdir = fileparts(path_prefix);
if ~isempty(pathdir) && ~exist(pathdir, 'dir')
  mkdir(pathdir);
end

% Give this matrix a name
p = length(A);
matname = ['net_', num2str(p), '_0X', BKDRHash(A)];
if issparse(A)
  matname = [matname '_sp'];
end
matpath = [path_prefix, matname, '.txt'];

% Solve file name collision by adding extra characters, if any
matname0 = matname;
k = 0;
while exist(matpath, 'file')
  stat = dir(matpath);
  if stat.bytes == 0
    break
  end
  B = load('-ascii', matpath);
  if issparse(A)
    B = spconvert(B);
  end
  if ~any(size(A) - size(B)) &&  isempty(find(A-B, 1))
    return
  end
  k = k + 1;
  matname = sprintf('%s_%s', matname0, lower(dec2base(k,16)));
  matpath = [path_prefix, matname, '.txt'];
  warning('save_network:hash', 'hash collision occured! File is renamed.');
end

if issparse(A)
  fd = fopen(matpath, 'w');
  if fd == -1
    error('Unable to open output file `%s''', matpath);
  end
  fprintf(fd, '# sparse network, size = %d x %d, format: i,j,value\n', size(A,1), size(A,2));
  [ii, jj, val] = find(A);
  ijv = [ii jj val];
  if all(floor(val) == val)
    % All data are integers
    fprintf(fd, '%d %d %d\n', ijv');
  else
    fprintf(fd, '%d %d %.17g\n', ijv');
  end
  if size(ijv,1) == 0
    fprintf(fd, '%d %d %d\n', [size(A),0]');
  end
  fclose(fd);
  return
end

% Save the matrix to the file
if all(floor(A(:)) == A(:))
  % Matrix is consist of integers, so let's save it in a cleaner way
  if max(A(:)) <= 9 && min(A(:)) >= 0 
    cA = char(zeros(size(A,1), 2*size(A,2)) + ' ');
    cA(:, 1:2:end) = A + '0';
    A = cA;
  else
    A = int2str(int32(A));      % most slow step
  end
  fd = fopen(matpath, 'w');
  if fd == -1
    error('Unable to open output file `%s''', matpath);
  end
  for jj=1:p
    fprintf(fd, '%s\n', A(jj,:));
  end
  fclose(fd);
else
  % Save the matrix in full precision ascii format
  save('-ascii', '-double', matpath, 'A');
end

%% test
%A = rand(10)<0.5;
%[matpath matname] = save_network(A, 'ab');
%B = get_network(matpath);
%assert(~any(A(:)-B(:)))

%[matpath matname] = save_network(A, 'ab/');
%B = get_network(matpath);
%assert(~any(A(:)-B(:)))

%A = sparse(A);
%[matpath matname] = save_network(A, './');
%B = get_network(matpath);
%assert(~any(A(:)-B(:)))

%tic; A=sparse(rand(1000)<0.5); toc; tic; fn = save_network(A, 'ab/'); toc
%B = get_network(fn);
%assert(~any(A(:)-B(:)))

%tic; A=sparse((rand(1000)<0.5).*rand(1000)); toc; tic; fn = save_network(A, 'ab/'); toc
%tic; B = get_network(fn); toc;
%assert(~any(A(:)-B(:)))

