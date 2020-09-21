function [] = writeDefFile(t0,tf,freetf,dims,disc_steps,disc_method,path)

fid = fopen(path,'w+');

% write input data
fprintf(fid,'# This file defines all dimensions and parameters\n# values for your problem :\n');

fprintf(fid,'\n# Initial and final time :\ntime.free string %s\n',freetf);
fprintf(fid,'time.initial double %s\ntime.final double %s\n',t0,tf);

fprintf(fid,'\n# Dimensions :\n');
labels={'boundarycond','state', 'control', 'algebraic', 'parameter', 'constraint', 'constant'};
for i=1:7
    fprintf(fid,'%s.dimension integer %d\n',labels{i},dims(i));
end

fprintf(fid,'\n# Discretization :\ndiscretization.steps integer %s\ndiscretization.method string %s\n\n'...
    ,disc_steps,disc_method);

% write fixed part
fixedpart = {'# Optimization :'
    'optimization.type string single'
    'batch.type integer 0'
    'batch.index integer 0'
    'batch.nrange integer 1'
    'batch.lowerbound double 0'
    'batch.upperbound double 0'
    'batch.directory string none'
    ''
    '# Initialization :'
    'initialization.type string from_init_file'
    'initialization.file string none'
    ''
    '# Parameter identification :'
    'paramid.type string false'
    'paramid.separator string ,'
    'paramid.file string no_directory'
    'paramid.dimension integer 0'
    ''
    '# Names :'
    ''
    '# Solution file :'
    'solution.file string problem.sol'
    ''
    '# Iteration output frequency :'
    'iteration.output.frequency integer 0'};

for line = fixedpart'
    fprintf(fid,'%s\n',line{1});
end
fclose(fid)




end

