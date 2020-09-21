function names = read_names(fid,dim)
names = cell(dim,1);
for i=1:dim
    names{i} = rs(fid);
end
end
