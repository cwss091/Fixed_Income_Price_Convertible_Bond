%% Calibrating Heston Model
% Copyright @ Yuting Wu in Georgia Tech
clc
clear all
close all

Put_Option = xlsread('Put_Option.xlsx',1,'C21:K26');
Time_to_Maturity = [12 37 57 122 187 452]./252;
Strike_Price = [20 30 35 40 42.5 45 47.5 50 55];
Put_Price = Put_Option';
r = 0.0114;

S_0 = 47.42;
v_0 = 0.214491667;

big = 10^10;
options = optimset('MaxFunEvals',big,'MaxIter',big); % I guess estimating nine values is burdensome. So, we need to increase the maximum number of function evaluations so that we can achieve the optimized value.

parameters_plugged = [r S_0 v_0];
f = @(x)Implied_Volatilty_Distance_Heston( x , parameters_plugged , Strike_Price , Time_to_Maturity , Put_Price);
initial = [0 0 0 -0.1];% starting value
[sol,fval,exitflag,output,grad,hessian] = fminunc(f,initial,options);

sigma_hat = 0.26313245 + 0.26313245*exp(sol(1))/(1+exp(sol(1)));
xi_hat = 0.885508883 + 0.885508883*exp(sol(2))/(1+exp(sol(2)));
kappa_hat = (xi_hat^2)/(2*(sigma_hat^2)) + exp(sol(3));
rho_hat = -0.1 - 0.9*exp(sol(4))/(1+exp(sol(4)));
