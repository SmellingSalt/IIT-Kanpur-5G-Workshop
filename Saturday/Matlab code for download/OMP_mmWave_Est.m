function h_b_omp = OMP_mmWave_Est(y,Qbar,thrld)
[rq,cq] = size(Qbar); %% obtain the number of rows and columns
set_I = zeros(cq,1);  %% will contain the active column index
r_prev = zeros(rq,1); %% initialize the previous residue
h_b_omp = zeros(cq,1);
r_curr = y; %% initialize the current residue
Q_a = []; 
ix1 = 1;
while(abs((norm(r_prev))^2 - (norm(r_curr))^2) > thrld)
    % obtain the index of the largest absolute inner product
    [m_val,m_ind] = max(abs(Qbar'*r_curr)); 
    set_I(ix1) = m_ind;
    Q_a = [Q_a, Qbar(:,m_ind)]; %% update the active column set
    H_b_ls = pinv(Q_a)*y; %% find the LS sol for active column set
    r_prev = r_curr;
    r_curr = y - Q_a*H_b_ls; %% find the approximation error
    ix1 = ix1 + 1;
end
% set the LS solution in the active indices
h_b_omp(set_I(1:ix1-1)) = H_b_ls;