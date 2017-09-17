% Get the adjacency matrix from its (file) name.
% Assume the matrix is stored in plain text file. Which essentially can be
%   loaded  by `network = load('-ascii', [netstr, '.txt'])'
% It will search current working dir (or path_prefix if specified) and
%   dir of this function.
% It also possible to specifies the full path in `netstr'.

function [network, matname] = get_network(netstr, path_prefix)
if ~exist('netstr','var')
  error('Usage: [network, matname] = get_network(netstr [, path_prefix])');
end

if ~ischar(netstr)
  error('Input should be the name of the matrix');
end

if ~exist('path_prefix','var')
  % Use default dir (usually working dir)
  path_prefix = '';
end
% always consider path_prefix as dir name
e = filesep;
if isempty(netstr)
  matname = '-';
  network = [1];
  return
end

pathdir0 = fileparts(mfilename('fullpath'));
% path candidates, priority from high to low
s_matname = {...
  [path_prefix, netstr, '.txt'],...
  [path_prefix, netstr],...
  [pathdir0, e, 'network', e, netstr, '.txt'],...
  [pathdir0, e, 'network', e, netstr]...
};
for matname = s_matname
  matname = matname{1};
  %fprintf('f=%s\n', matname);
  if exist(matname, 'file')
    network = load('-ascii', matname);
    break
  end
end
if ~exist('network','var')
  warning('Network file not found!');
  network = [];
  matname = [];
end

if size(network, 2) == 3 && size(network, 1) != 3
  % This is a sparse matrix
  network = spconvert(network);
end

end
