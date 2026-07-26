% Neoclassical Growth Model with Capital and Exogenous Labor Supply (Nonlinear version)
close all;

var

C // consumption
K // capital
R // rental cost of capital
W // wage
Z // TFP: This Z is actually exp(z) in the handout.
;

parameters 
alpha beta delta rho;

alpha = 0.34; // Share of capital income 
beta = 0.96; // Time discount factor
delta = 0.04; // Depreciation rate
rho = 0.9; // Persistence of TFP shock

varexo
eps_z;


shocks;
var eps_z; stderr 0.5;
end;

model;

1 / C = (beta * (R(+1) + 1 - delta)) / C(+1);
R = alpha * Z * K(-1)^(alpha - 1);
W = (1-alpha) * Z * K(-1)^(alpha);
C + K = (R + 1 - delta) * K(-1) + W;
log(Z) = rho * log(Z(-1)) + eps_z;

end;

steady_state_model;

Z = 1;
R = 1 / beta - 1 + delta;
K = (1 / alpha)^(1/(alpha-1)) * (1 / beta - 1 + delta)^(1 / (alpha - 1));
C = K^alpha - delta * K;
W = (1-alpha) * K^alpha;

end;

steady; // Dynare can also solve for the steady state!
resid;
check;
model_diagnostics;

stoch_simul(order=1, irf=20, loglinear) K C W;
stoch_simul(order=1, irf=20) R;
