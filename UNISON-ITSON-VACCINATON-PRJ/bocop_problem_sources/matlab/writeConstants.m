function [] = writeConstants(constants,path)

dimConstants = size(constants,1);
fid = fopen(path,'w+');
fprintf(fid,'# This file contains the values of the constants of your problem.\n');
fprintf(fid,'# Number of constants used in your problem :\n%d\n\n',dimConstants);
fprintf(fid,'# Values of the constants : \n');
fprintf(fid,'%f\n',constants);
fclose(fid);

end

