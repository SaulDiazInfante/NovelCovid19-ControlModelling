% load parameters and initial conditions
fname = 'vaccine_multidoses_parameters.json';
fid = fopen(fname);
raw = fread(fid, inf);
str = char(raw'); 
fclose(fid);
prm = jsondecode(str);
u_zero = [prm.S_0, prm.E_0, prm.I_S_0, prm.I_A_0, ...
            prm.H_0, prm.R_0, prm.D_0, ...
            prm.X_S_0, prm.X_E_0, prm.X_I_A_0, prm.X_R, ...
            prm.Y_S_0, prm.Y_E_0, prm.Y_I_A_0, prm.Y_R, ...         
            prm.X_VAC_1_A_0, prm.X_VAC_2_A_0];

opts=odeset('reltol',1e-8, 'maxstep', 0.1);
[t, sol] = ode45(@vaccine_multidoses_rhs, ...
    0:1:prm.T, u_zero, opts, prm);         

s = sol(:, 1);
e = sol(:, 2);
i_s = sol(:, 3);
i_a = sol(:, 4);
h = sol(:, 5);
r = sol(:, 6);
d = sol(:, 7);
% First dose vaccine type A
x_s =  sol(:, 8);
x_e = sol(:, 9);
x_i_a = sol(:, 10);
x_r= sol(:, 11);
% Second dose vaccine type B
y_s = sol(:, 12);
y_e = sol(:, 13);
y_i_a = sol(:, 14);
y_r = sol(:, 15);
x_vac_1_a = sol(:, 16);
x_vac_2_a = sol(:, 17);
% Merged compartments
ss = s + x_s + y_s;
ee = e + x_e + y_e;
ii_a = i_a + x_i_a + y_i_a;
rr = r + x_r + y_r;
%
vac_1_doses_t = prm.lambda_v_1_a * (s + e + i_a + r);
vac_2_doses_t = prm.lambda_v_2_a * (y_s + y_e + y_i_a + y_r);
%
tiledlayout('flow')
nexttile
plot(t, ss, 'DisplayName', 'Whole S')
hold on
plot(t, s, 'DisplayName', 'S')
plot(t, x_s, 'DisplayName', 'X_s' )
plot(t, y_s, 'DisplayName', 'Y_s' )
ylabel('Suceptible') 
legend
%
nexttile
plot(t, e)
hold on
plot(t, x_e, 'DisplayName', 'X_E' )
plot(t, y_e, 'DisplayName', 'Y_E' )
ylabel('Exposed') 
legend
%
nexttile
plot(t, i_s)
ylabel('Symptomatic') 
%
nexttile
plot(t, i_a, 'DisplayName', 'I_A' )
hold on
plot(t, x_i_a, 'DisplayName', 'X_I_A' )
plot(t, y_i_a, 'DisplayName', 'Y_I_S' )
ylabel('Asymptomatic') 
legend
%
nexttile
plot(t, h)
ylabel('Hospitalized') 
%
nexttile
plot(t, rr, 'DisplayName', 'Whole R' )
hold on
plot(t, r, 'DisplayName', 'R' )
plot(t, x_r, 'DisplayName', 'X_R' )
plot(t, y_r, 'DisplayName', 'Y_R' )
ylabel('Recover') 
%
nexttile
plot(t, d)
ylabel('Death') 
%
nexttile
plot(t, x_vac_1_a)
hold on
plot(t, x_vac_2_a, '--r')
ylabel('Doses') 
xlim([0, 100])
%
nexttile
plot(t, x_vac_1_a)
hold on
plot(t, vac_1_doses_t)
ylabel('First Dose')
xlim([0, 365])
%
nexttile
plot(t,  x_vac_2_a)
hold on
plot(t, vac_2_doses_t)
ylabel('Second Dose')
xlim([0, 365])
