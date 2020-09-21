function [] = writeIpopt(options_ipopt,path)

fid = fopen(path,'w+');
for line = options_ipopt'
   fprintf(fid,'%s %s\n',line{1},line{2}); 
end
fclose(fid);    

end

