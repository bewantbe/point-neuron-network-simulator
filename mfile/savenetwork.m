% Save the adjacency matrix `A' in path `pathdir' with a hashed file name:
%   matpath = savenetwork(A, pathdir);
% To read the matrix from the file:
%   A = load('-ascii', matpath);

function [matpath, matname] = savenetwork(A, pathdir)

% sane test of the matrix A
if ~exist('A','var') || isempty(A) || (~isnumeric(A) && ~islogical(A))
  error('Input a numeric matrix please!');
end
if ~ismatrix(A) || size(A,1) ~= size(A,2)
  error('Not a square matrix!');
end
A = double(A);  % int32,logical to double, so get consistent hash result

e = filesep;
if ~exist('pathdir','var')
  pathdir = '';    % save to default dir (usually working dir)
end
% Always consider `pathdir' as dir name
if ~isempty(pathdir) && pathdir(end) ~= '/' && pathdir(end) ~= e
  pathdir = [pathdir e];
end
if ~exist(fileparts(pathdir), 'dir')
  mkdir(pathdir);
end

% Give this matrix a name
p = length(A);
matname = ['net_', num2str(p), '_0X', BKDRHash(mat2str(A))];
matpath = [pathdir, matname, '.txt'];

% Solve file name collision by adding extra characters, if any
matname0 = matname;
k = 0;
while exist(matpath, 'file')
  B = load('-ascii', matpath);
  if norm(A(:)-B(:)) == 0
    return
  end
  k = k + 1;
  matname = sprintf('%s_%s', matname0, lower(dec2base(k,16)));
  matpath = [pathdir, matname, '.txt'];
  warning('savenetwork:hash', 'hash collision occured! File is renamed.');
end

% Save the matrix to the file
if all(floor(A(:)) == A(:))
  % Matrix is consist of integers, so let's save it in a cleaner way
  A = int2str(int32(A));
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
