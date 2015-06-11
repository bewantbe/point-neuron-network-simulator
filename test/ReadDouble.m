%

function X = ReadDouble(fname, p)
if p<=0
  error('');
end
fid = fopen(fname);
X = fread(fid, [p Inf], 'double');
fclose(fid);
