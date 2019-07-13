function [BB, RF] = SOMP_mmW_precoder(Opt, Dict, Ryy, numRF)
[rq, cq] = size(Dict); % size of dictionary
Res = Opt; % initial Residue
RF = [];   % RF precoder/ combiner columns
for iter1 = 1:numRF
    phi = Dict'*Ryy*Res; % finding projection
    phi_phiH = phi*phi';
    for iter2 = 1:cq;
        phi_phiH_diag(iter2) = phi_phiH(iter2,iter2);
    end
    % choosing maximum projection
    [MaxVal,MaxInd] = max(phi_phiH_diag);    
    % augmenting RF precoder/ combiner matrix
    RF = [RF Dict(:,MaxInd)];  
    % calculate BB precoder/ combiner
    BB = inv(RF'*Ryy*RF)*RF'*Ryy*Opt;
    % updating the residue
    Res = (Opt-RF*BB)/(norm((Opt-RF*BB),'fro'));    
end
