%Tail Probability Curves
clear; clc; rng(2026);

%Parameters
S0 = 100;
mu = 0.08;

T  = 1.0;
N  = 252;
dt = T / N;

M = 30000;   % number of Monte Carlo paths

%GBM
sigma_gbm = 0.20;

%Heston
sigma_base = 0.20;
v0    = sigma_base^2;
kappa = 2.0;
theta = sigma_base^2;
xi    = 0.40;
rho   = -0.50;

%Merton
sigma_m = 0.20;
lambda = 1.0;
muJ    = -0.10;
sigmaJ = 0.30;

% Drift for Merton
EJ = exp(muJ + 0.5*sigmaJ^2);
mu_merton = mu - lambda*(EJ - 1);

%Storee
SgT = zeros(M,1);
ShT = zeros(M,1);
SmT = zeros(M,1);

%Simulations
for i = 1:M
    Sg = S0;
    Sh = S0; v = v0;
    Sm = S0;

    for n = 1:N
       
        Z1 = randn;
        Z2 = randn;

        %GBM
        Sg = Sg * exp((mu - 0.5*sigma_gbm^2)*dt ...
              + sigma_gbm*sqrt(dt)*Z1);

        %Heston
        Z2corr = rho*Z1 + sqrt(1 - rho^2)*Z2;
        v = v + kappa*(theta - v)*dt ...
              + xi*sqrt(max(v,0))*sqrt(dt)*Z2corr;
        v = max(v,0);
        Sh = Sh * exp((mu - 0.5*v)*dt ...
              + sqrt(v)*sqrt(dt)*Z1);

        %Merton
        Zm = randn;
        Sm_cont = Sm * exp((mu_merton - 0.5*sigma_m^2)*dt ...
                  + sigma_m*sqrt(dt)*Zm);

        K = poissrnd(lambda*dt);
        if K > 0
            J = exp(muJ + sigmaJ*randn(K,1));
            Sm = Sm_cont * prod(J);
        else
            Sm = Sm_cont;
        end
    end

    SgT(i) = Sg;
    ShT(i) = Sh;
    SmT(i) = Sm;
end

%Log-returns
R_gbm = log(SgT / S0);
R_hes = log(ShT / S0);
R_mer = log(SmT / S0);

%Tail prob curves
xVals = linspace(0, 0.6, 121); 
tail_gbm = zeros(size(xVals));
tail_hes = zeros(size(xVals));
tail_mer = zeros(size(xVals));

absR_gbm = abs(R_gbm);
absR_hes = abs(R_hes);
absR_mer = abs(R_mer);

for k = 1:length(xVals)
    x = xVals(k);
    tail_gbm(k) = mean(absR_gbm > x);
    tail_hes(k) = mean(absR_hes > x);
    tail_mer(k) = mean(absR_mer > x);
end


figure;
plot(xVals, tail_gbm, 'k--', 'LineWidth', 2); hold on;
plot(xVals, tail_hes, 'r-',  'LineWidth', 2);
plot(xVals, tail_mer, 'b-',  'LineWidth', 2);

xlabel('x', 'FontSize', 16);
ylabel('P(|R| > x)', 'FontSize', 16);
title('Tail Probability Curve for Terminal Log-Returns', 'FontSize', 18);

legend({'GBM','Heston','Merton'}, 'FontSize', 28, 'Location', 'northeast');
set(gca, 'FontSize', 28);
grid on; box on;
set(gcf, 'Color', 'w');
