%Merton Price Paaths
clear; clc; rng(2026);

S0 = 100;
mu = 0.08;

T  = 1.0;
N  = 252;
dt = T / N;
time = linspace(0, T, N+1);

nPlot = 10;   % number of paths

sigma  = 0.20;     % diffusion volatility
lambda = 1.0;      % expected jumps per year
muJ    = -0.10;    % mean log jump size
sigmaJ = 0.30;     % standard log jump size

EJ = exp(muJ + 0.5*sigmaJ^2);     
mu_adj = mu - lambda*(EJ - 1);

S_mer = zeros(nPlot, N+1);
S_mer(:,1) = S0;

for i = 1:nPlot
    S = S0;

    for n = 1:N
        Z = randn;
        S_cont = S * exp((mu_adj - 0.5*sigma^2)*dt ...
                 + sigma*sqrt(dt)*Z);

        % Jump
        K = poissrnd(lambda*dt);
        if K > 0
            J = exp(muJ + sigmaJ*randn(K,1));
            S = S_cont * prod(J);
        else
            S = S_cont;
        end

        S_mer(i, n+1) = S;
    end
end

figure;
plot(time, S_mer.', 'LineWidth', 1.3);

xlabel('Time (years)', 'FontSize', 16);
ylabel('S_t', 'FontSize', 16);
title('Sample Price Paths: Merton Jump-Diffusion Model', 'FontSize', 28);

set(gca, 'FontSize', 26);
box on;
set(gcf, 'Color', 'w');
