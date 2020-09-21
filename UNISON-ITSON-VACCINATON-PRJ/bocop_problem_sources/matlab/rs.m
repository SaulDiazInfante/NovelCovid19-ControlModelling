function result = rs(fid)
    result = fscanf(fid,'# %*s %*s %s \n',1);
end