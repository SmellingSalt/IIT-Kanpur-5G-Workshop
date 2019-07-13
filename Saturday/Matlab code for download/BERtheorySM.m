% Analytical Expression [Jeganathan et al, "Spatial Modulation: Optimal 
% Detection and Performance Analysis", IEEE Comm. Letters, Aug 2008]
function BERtheory = BERtheorySM(SNRdB, Nt, Nr, M)
    rho = 10^(SNRdB/10);
    % sigma-square alpha
    sigma2alpha = rho/2;
    % muAlpha
    muAlpha = 0.5*(1 - sqrt(sigma2alpha/(1+sigma2alpha)));
    anSum = 0;
    for K = 0:Nr-1
        anSum = anSum + nchoosek(Nr-1+K,K)*(1-muAlpha)^K;
    end
    BERtheory = 2*Nt*(muAlpha^Nr)*anSum/M;
end