%Steffen Docken
%27-7-22
%Code to calculate AIC for a generic model given:
% lnL: maximum log-likelihood 
% k: number of parameters
% n: number of measurements
% Corrected: 0 (uncorrected AIC), 1 (corrected AIC)

function AIC_val = AIC_calc_generic(lnL, k, n, Corrected)

AIC_uncorrect = 2*k - 2*lnL;

if (Corrected == 0)
    AIC_val = AIC_uncorrect;
else
    AIC_val = AIC_uncorrect + (2*k^2 +2*k)/(n - k - 1);
end