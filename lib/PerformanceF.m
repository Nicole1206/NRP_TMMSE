function [SINR,SOI] = PerformanceF(EnerS,Q,R,R0,alpha,ad,adS)

Wopt = (inv(R)^alpha * ad) / (ad' * inv(R)^alpha * ad);


% Ws = 10 * log10(abs(Wopt' * Rss * Wopt));
% Wgi = 10 * log10(abs(Wopt' * Rii1 * Wopt) + abs(Wopt' * Rii2 * Wopt) + abs(Wopt' * Rww * Wopt));
SOI = 10 * log10(abs(Wopt' * R0 * Wopt));
% SINR = Ws - Wgi;
SINR = 10 * log10(10^(EnerS/10) * abs(Wopt' * adS)^2 / real(Wopt' * Q * Wopt));