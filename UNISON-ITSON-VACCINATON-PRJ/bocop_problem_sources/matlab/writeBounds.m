function [] = writeBounds(dims,bounds,path)

fid = fopen(path,'w+');
fprintf(fid,'# This file contains the bounds of your problem.\n');
fprintf(fid,'# Bounds are stored in standard format :\n# [lower bound]  [upper bound] [type of bound]\n\n');
fprintf(fid,'# Dimensions (i&f conditions, y, u, z, p, path constraints) :\n');
fprintf(fid,'%d ',dims(1:end-1));
fprintf(fid,'\n\n');

dim_bounds = size(bounds,1)
for i=1:dim_bounds
    lb = bounds(i,1);
    ub = bounds(i,2);
    if (lb == ub)
        type = 'equal';
    elseif (lb == -2e20 && ub == 2e20)
        type = 'free';
    elseif (lb > -2e20 && ub < 2e20)
        type = 'both';
    elseif (lb > -2e20)
        type = 'lower';
    else
        type = 'upper';
    end
    fprintf(fid,'%f %f %s\n',lb,ub,type);
end

fclose(fid);

end

