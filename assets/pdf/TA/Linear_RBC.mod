% Neoclassical Growth Model with Capital and Exogenous Labor Supply (Linear version)
close all;

var

C // consumption
K // capital
R // rental cost of capital
W // wage
z // TFP
;

parameters 
alpha beta delta rho barK barR barC barW;

alpha = 0.34; // Share of capital income 
beta = 0.96; // Time discount factor
delta = 0.04; // Depreciation rate
rho = 0.9; // Persistence of TFP shock
barK = (1 / alpha)^(1/(alpha-1)) * (1 / beta - 1 + delta)^(1 / (alpha - 1));
barR = 1 / beta - 1 + delta;
barC = barK^alpha - delta * barK;
barW = (1-alpha) * barK^alpha;

varexo
eps_z;

shocks;
var eps_z; stderr 0.5;
end;

model;

- C = - C(+1) + beta * barR * R(+1);
R = z + (alpha - 1) * K(-1);
W = z + alpha * K(-1);
barC * C + barK * K = barK * barR * R + (barR + 1 - delta) * barK * K(-1) + barW * W;
z = rho * z(-1) + eps_z;

end;

steady_state_model;

z = 0;
C = 0;
R = 0;
K = 0;
W = 0;

end;

steady; // Dynare can also solve for the steady state!
resid;
check;
model_diagnostics;

stoch_simul(order=1, irf=20) K C R W;
