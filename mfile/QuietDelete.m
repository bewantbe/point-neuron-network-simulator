function QuietDelete(fname)
if exist(fname, 'file')
  delete(fname);
end
