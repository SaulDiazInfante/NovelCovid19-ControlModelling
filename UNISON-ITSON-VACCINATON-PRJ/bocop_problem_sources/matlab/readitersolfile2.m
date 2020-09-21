% Matlab script to read bocop .sol.iter file
% Pierre Martinon
% Inria and CMAP Ecole Polytechnique
% 2015-2017


% TODO: preallocate vectors to speedup read

function [dim_state,dim_control,dim_algebraic,dim_optimvars,dim_constant,...
    dim_boundarycond,dim_pathcond,time_steps,...
    time,stage,state,control,algebraic,constants,optimvars,...
    state_lb,state_ub,control_lb,control_ub,...
    objective,constraints_Infnorm] = readitersolfile2(path,verbose)


%open file
disp('');
[fid] = fopen(path,'r');
if fid==-1
   fprintf('Incorrect path for solution file: %s\n', path);
   return
elseif (verbose > 0)
   fprintf('Reading file ... %s\n',path);
end

%%%%%%%%%%%%%%%%%%%%%
% read .def file copy
for i=1:14
    blob = fgetl(fid);
end

% Initial and final time block
free_time = rs(fid);
initial_time = rf(fid);
final_time = rf(fid);

for i=1:2
    blob = fgetl(fid);
end

% Dimension block
dim_state = ri(fid);
dim_control = ri(fid);
dim_algebraic = ri(fid);
dim_optimvars = ri(fid);
dim_constant = ri(fid);
dim_boundarycond = ri(fid);
dim_pathcond = ri(fid);
if (verbose > 0)
fprintf('Dimensions state: %d control: %d algebraic: %d parameter: %d\n'...
    , dim_state,dim_control,dim_algebraic,dim_optimvars);
fprintf('Dimensions constant: %d boundarycond: %d constraint: %d\n' ...
    , dim_constant, dim_boundarycond, dim_pathcond);
end

for i=1:2
    blob = fgetl(fid);
end

% Discretization block
time_steps = ri(fid);
disc_method = rs(fid);
if (verbose > 0)
    fprintf('Discretization: %s    Time steps: %d\n',disc_method,time_steps);
end

for i=1:2
    blob = fgetl(fid);
end

% Optimization block
optimization_type = rs(fid);
blob = fgetl(fid);
%test for old files without batch type !
if (strfind(blob,'type'))
batch_type =  sscanf(blob,'# batch.type integer %d');
batch_index = ri(fid);
else
batch_index = sscanf(blob,'# batch.index integer %d');
end
batch_nrange = ri(fid);
batch_lowerbound = rf(fid);
batch_upperbound = rf(fid);
batch_directory = rs(fid);

for i=1:2
    blob = fgetl(fid);
end

% Initialization block
initialization_type = rs(fid);
initialization_file = rs(fid);

for i=1:2
    blob = fgetl(fid);
end

% Parameter identification block
%test for old files without parameter identification block !
if (strfind(blob,'Parameter'))
paramid_type = rs(fid);
paramid_separator = rs(fid);
paramid_file = rs(fid);
blob = fgetl(fid);
if (strfind(blob,'paramid'))
paramid_dimension =  sscanf(blob,'# paramid.dimension integer %d');
blob = fgetl(fid);
end
blob = fgetl(fid);
end

% Names block
state_names = read_names(fid,dim_state);
control_names = read_names(fid,dim_control);
algebraic_names = read_names(fid,dim_algebraic);
parameters_names = read_names(fid,dim_optimvars);
boundarycond_names = read_names(fid,dim_boundarycond);
constraint_names = read_names(fid,dim_pathcond);
constant_names = read_names(fid,dim_constant);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read .bounds
line = fgetl(fid);
while ~strcmp('# ** problem.bounds',deblank(line))
    line = fgetl(fid);
end
for i=1:10
    blob = fgetl(fid);
end

%boundary conditions bounds
boundarycond_lb=[];
boundarycond_ub=[];
for i=1:dim_boundarycond
   boundarycond_lb(i) = fscanf(fid,'# %f',1);
   boundarycond_ub(i) = fscanf(fid,'%f',1);
   line = fgetl(fid);
end
for i=1:2
    blob = fgetl(fid);
end

%state bounds
state_lb=[];
state_ub=[];
for i=1:dim_state
   state_lb(i) = fscanf(fid,'# %f',1);
   state_ub(i) = fscanf(fid,'%f',1);
   line = fgetl(fid);
end
for i=1:2
    blob = fgetl(fid);
end

%control bounds
control_lb=[];
control_ub=[];
for i=1:dim_control
   control_lb(i) = fscanf(fid,'# %f',1);
   control_ub(i) = fscanf(fid,'%f',1);
   line = fgetl(fid);
end
for i=1:2
    blob = fgetl(fid);
end

%algebraic
algebraic_lb=[];
algebraic_ub=[];
for i=1:dim_algebraic
   algebraic_lb(i) = fscanf(fid,'# %f',1);
   algebraic_ub(i) = fscanf(fid,'%f',1);
   line = fgetl(fid);
end
for i=1:2
    blob = fgetl(fid);
end

%optimvars
optimvars_lb=[];
optimvars_ub=[];
for i=1:dim_optimvars
   optimvars_lb(i) = fscanf(fid,'# %f',1);
   optimvars_ub(i) = fscanf(fid,'%f',1);
   line = fgetl(fid);
end
for i=1:2
    blob = fgetl(fid);
end

%pathcond bounds
pathcond_lb=[];
pathcond_ub=[];
for i=1:dim_pathcond
   pathcond_lb(i) = fscanf(fid,'# %f',1);
   pathcond_ub(i) = fscanf(fid,'%f',1);
   line = fgetl(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read .constants
line = fgetl(fid);
while ~strcmp('# ** problem.constants',deblank(line))
    line = fgetl(fid);
end

for i=1:7
    blob = fgetl(fid);
end
constants = zeros(1,dim_constant);
for i=1:dim_constant
constants(i) = fscanf(fid,'# %f \n',1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%skip initialization variables 

%get real number of steps
line = fgetl(fid);
while size(strfind(deblank(line),'merge')) == 0  &  ~strcmp('# *****     SOLUTION     *****',deblank(line))  
    line = fgetl(fid);
end

if (strfind(line,'discretization.steps.after.merge'))
time_steps2 = sscanf(line,'# discretization.steps.after.merge %d');
if (time_steps2 ~= time_steps)
    if (verbose > 0)
        fprintf('Actual number of time steps after merge: %d\n',time_steps2);
    end
    time_steps = time_steps2;
end
line = fgetl(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read solution
for i=1:6
    blob = fgetl(fid);
end

tline = fgets(fid);
objective=fscanf(fid,'%f',1);
if (verbose > 0)
    fprintf('Objective value: %f\n',objective)
end
tline = fgets(fid);
tline = fgets(fid);
constraints_L2norm=fscanf(fid,'%f',1);
tline = fgets(fid);
tline = fgets(fid);
constraints_Infnorm=fscanf(fid,'%f',1);
tline = fgets(fid);
tline = fgets(fid);
time_stage=fscanf(fid,'%f',1);

%times
alltimes=[]; borne=[]; ntime=[];
alltimes = fscanf(fid,'%f',(time_steps+1) + time_steps*time_stage); 
borne=length(alltimes);
ntime = alltimes(1:1+time_stage:end);
p=0;
nstage=[];
for i=1:(time_stage+1):(borne-time_stage)
    for n=1:time_stage, p=p+1; indice=(i+n); nstage(p)=alltimes(indice); end
end
nstage = nstage';
for i=1:3, tline = fgets(fid); end

%state
state=[];
for i=1:dim_state
 state(1:time_steps+1,i) = fscanf(fid,'%f',time_steps+1);
 for j=1:3, tline = fgets(fid); end
end

%control
control=[];
for i=1:dim_control
 control(1:time_steps*time_stage,i) = fscanf(fid,'%f',time_steps*time_stage);
 for j=1:3, tline = fgets(fid); end
end

%algebraic
algebraic=[];
for i=1:dim_algebraic
 algebraic(1:time_steps*time_stage,i) = fscanf(fid,'%f',time_steps*time_stage);
 for j=1:3, tline = fgets(fid); end
end

%parameter
optimvars=[];
if dim_optimvars~=0
 optimvars(1:dim_optimvars) = fscanf(fid,'%f',dim_optimvars);
 tline = fgets(fid); 
end
for i=1:2, tline = fgets(fid); end

% %boundary conditions
% defboundary=[];
% boundarycond=[];
%  for i=1:dim_boundarycond
%  tline = fgets(fid);
%  defboundary{i} = textscan(tline,'%s %s %s','delimiter',' ');
%  boundarycond(i)=str2double(defboundary{i}{2}{1});
%  end
%  for i=1:2, tline = fgets(fid); end;
% 
% %path_constraint
% path_constraint=[];
% for i=1:dim_pathcond
%  tline=fgets(fid);
%  path_constraint(1:time_steps*time_stage,i) = fscanf(fid,'%f',time_steps*time_stage);
%  for i=1:3, tline = fgets(fid); end
% end
% 
% %dynamic_constraint
% dynamic_constraint=[];
%  for i=1:dim_state
%   tline=fgets(fid);
%   dynamic_constraint(1:time_steps,i) = fscanf(fid,'%f',time_steps);
%   for j=1:3, tline = fgets(fid); end
%  end
%  dynamic_constraint;

% %skipping boundary multipliers comments
% tline = fgets(fid);
% for i=1:4, tline = fgets(fid); end

% %boundarycond_mult
% boundarycond_mult=[];
% if dim_boundarycond~=0
%  boundarycond_mult(1:dim_boundarycond)=fscanf(fid,'%f',dim_boundarycond);
%  tline = fgets(fid);
% end
% for i=1:2, tline = fgets(fid); end 

% %pathconstraint_mult
% pathconstraint_mult=[];
% for i=1:dim_pathcond
%  pathconstraint_mult(1:time_steps*time_stage,i) = fscanf(fid,'%f',time_steps*time_stage);
%  for j=1:3, tline = fgets(fid); end
% end
% 
% %adjoint_state (dynamic multiplier)
% adjoint=[];
% for i=1:dim_state
%  adjoint(1:time_steps,i) = fscanf(fid,'%f',time_steps);
%  for j=1:3, tline = fgets(fid); end
% end

% %coefficients of discretization method
% for i=1:2, tline = fgets(fid); end
% coef_a=[];
% coef_a(1:time_stage*time_stage)=fscanf(fid,'%f',time_stage*time_stage);
% for i=1:3, tline = fgets(fid); end
% coef_b=[];
% coef_b(1:time_stage)=fscanf(fid,'%f',time_stage);
% for i=1:3, tline = fgets(fid); end
% coef_c=[];
% coef_c(1:time_stage)=fscanf(fid,'%f',time_stage);
% for i=1:5, tline = fgets(fid); end

% %zl_states
% zl_state=[];
% for i=1:dim_state
%  zl_state(1:time_steps+1,i) = fscanf(fid,'%f',time_steps+1);
%  for j=1:3, tline = fgets(fid); end
% end
% 
% %zl_control
% zl_control=[];
% for i=1:dim_control
%  zl_control(1:time_steps*time_stage,i) = fscanf(fid,'%f',time_steps*time_stage);
%  for j=1:3, tline = fgets(fid); end
% end
% 
% %zl_algebraic
% zl_algebraic=[];
% for i=1:dim_algebraic
%  zl_algebraic(1:time_steps*time_stage,i) = fscanf(fid,'%f',time_steps*time_stage);
%  for j=1:3, tline = fgets(fid); end
% end
% 
% %zl_parameter
% zl_parameter=[];
% if dim_optimvars~=0
%  zl_parameter(1:dim_optimvars) = fscanf(fid,'%f',dim_optimvars);
%  tline = fgets(fid); 
% end
% for i=1:2, tline = fgets(fid); end  
% 
% 
% %zu_states
% zu_state=[];
% for i=1:dim_state
%  zu_state(1:time_steps+1,i) = fscanf(fid,'%f',time_steps+1);
%  for j=1:3, tline = fgets(fid); end
% end
% 
% %zu_control
% zu_control=[];
% for i=1:dim_control
%  zu_control(1:time_steps*time_stage,i) = fscanf(fid,'%f',time_steps*time_stage);
%  for j=1:3, tline = fgets(fid); end
% end
% 
% %zu_algebraic
% zu_algebraic=[];
% for i=1:dim_algebraic
%  zu_algebraic(1:time_steps*time_stage,i) = fscanf(fid,'%f',time_steps*time_stage);
%  for j=1:3, tline = fgets(fid); end
% end
% 
% %zu_parameter
% zu_parameter=[];
% zu_parameter(1:dim_optimvars) = fscanf(fid,'%f',dim_optimvars);
% for i=1:3, tline = fgets(fid); end
% 
% %cpu time
% cpu_time = fscanf(fid,'%f',1);
% for i=1:3, tline = fgets(fid); end
% 
% %iterations
% iterations = fscanf(fid,'%f',1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computed values

% %control_average
% control_average=[];
% for i=1:dim_control
% 	for j=1:time_steps
% 	    control_average(j,i) = 0; ind = 0;
% 	    for k=1:time_stage, ind = (j-1)*time_stage+k;
%             control_average(j,i) = control_average(j,i) + control(ind,i)*coef_b(k); 
%         end
%     end
% end

%non normalized time
if strcmp(free_time,'final')
time=initial_time+(optimvars(end)-initial_time).*ntime; 
else, time=initial_time+(final_time-initial_time).*ntime;
end

%non normalized stage
if strcmp(free_time,'final')
stage=initial_time+(optimvars(end)-initial_time).*nstage; 
else, stage=initial_time+(final_time-initial_time).*nstage;
end

