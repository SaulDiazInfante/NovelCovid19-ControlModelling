% load parameters and initial conditions
fname = 'vaccine_multidoses_parameters.json';
fid = fopen(fname);
raw = fread(fid, inf);
str = char(raw'); 
fclose(fid);
prm = jsondecode(str);
u_zero = [prm.S_0, prm.E_0, prm.I_S_0, prm.I_A_0, ...
            prm.H_0, prm.R_S_0, prm.R_A_0, prm.D_0, ...
            prm.X_S_0, prm.X_E_0, prm.X_I_A_0, prm.X_R_A, ...
            prm.Y_S_0, prm.Y_E_0, prm.Y_I_A_0, prm.Y_R_A, ...
            prm.X_VAC_1_A_0, prm.X_VAC_2_A_0];

opts=odeset('reltol',1e-8, 'maxstep', 0.1);
[t, sol] = ode45(@vaccine_multidoses_rhs, ...
    0:1:prm.T, u_zero, opts, prm);         

s = sol(:, 1);
e = sol(:, 2);
i_s = sol(:, 3);
i_a = sol(:, 4);
h = sol(:, 5);
r_s = sol(:, 6);
r_a = sol(:, 7);
d = sol(:, 8);
% First dose vaccine type A
x_s =  sol(:, 9);
x_e = sol(:, 10);
x_i_a = sol(:, 11);
x_r_a= sol(:, 12);
% Second dose vaccine type B
y_s = sol(:, 13);
y_e = sol(:, 14);
y_i_a = sol(:, 15);
y_r_a = sol(:, 16);
x_vac_1_a = sol(:, 17);
x_vac_2_a = sol(:, 18);

ss = s + x_s + y_s;
ee = e + x_e + y_e;
ii_a = i_a + x_i_a + y_i_a;
rr_a = x_r_a + y_r_a;

plot(t, i_s)

