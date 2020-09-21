function [] = writeInit(label,type,value,path)

dimValue = size(value);
fid = fopen([path label '.init'],'w+');

if (strcmp(label,'optimvars'))
    fprintf(fid,'# Optimization parameters starting point file.\n');
    fprintf(fid,'# This file contains initialization values\n');
    fprintf(fid,'# for all optimization parameters\n\n');
    fprintf(fid,'# Number of optimization parameters :\n%d\n\n',dimValue);
    fprintf(fid,'# Initial values : \n');
    fprintf(fid,'%f\n',value);
    fclose(fid);
    
else
    fprintf(fid,'# Starting point file.\n');
    fprintf(fid,'# This file contains initialization values\n');
    fprintf(fid,'# for variable %s\n\n',label);
    fprintf(fid,'# Type of initialization :\n%s\n\n',type);
    fprintf(fid,'# Constant value for the starting point : \n');
    fprintf(fid,'%f\n',value);
    fclose(fid);
    
end

end

