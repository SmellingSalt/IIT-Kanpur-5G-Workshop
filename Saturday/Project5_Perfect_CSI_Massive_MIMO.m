clear;
close all;
%% INITIALISATION
ITER = 1000;
K=10; % number of users
Mv =20:30:500; % number of BS antennas
Eu_dB= 10; Eu=10^(Eu_dB/10); %Total Power constraint in dB
rate_MRC=zeros(1,length(Mv));%This holds the rate of the MRC channel
bound_MRC=zeros(1,length(Mv));%Theoretical bound for MRC taken from some paper
rate_ZF = zeros(1,length(Mv));%This holds the rate of the ZF channel
bound_ZF =zeros(1,length(Mv));%Theoretical bound for ZF taken from some paper
%SINR MRC (NOTES)
%SINR ZF (NOTES)
%% SIMULATION
for it=1:ITER
    %D is the large scale fading coefficient diagonal
    %matrix(Beta1,...,BetaL for K users)
    D = Dmatrix(K); beta=diag(D); 
    for mx=1:length(Mv)
        %% Current base station under consideration
        M =Mv(mx);
        %% POWER SCALING
        %pu = Eu; % no power scaling
        pu = Eu/M;% power scaling
        %% CREATING THE CHANNEL AT THE RECIEVER
        H = sqrt(1/2)*(randn(M,K)+1i*randn(M,K));
        G =H*sqrt(D);
        %% SIMULATING
        for k= 1:K
            gk =G(:,k); % column of user k
            %% CONSTRUCTING THE MAXIMAL RATIO COMBINER FOR THE kth USER
            % NORM(gk)^4/NORM
            nr_MRC = pu*norm(gk)^4;
            dr_MRC = norm(gk)^2;
            nr_bound_MRC =pu*(M-1)*beta(k);%Formula for bound (PAPER)
            dr_bound_MRC = 1;
            for iu=1:K
                if (iu~=k) %(To compute the interfering user's power)
                    dr_MRC= dr_MRC+pu*abs(gk'*G(:,iu))^2;
                    dr_bound_MRC=dr_bound_MRC+pu*beta(iu);
                end
            end
            %% RATE FOR MRC
            rate_MRC(mx)=rate_MRC(mx)+log2(1+nr_MRC/dr_MRC);
            %% BOUND FOR MRC
            bound_MRC(mx)=bound_MRC(mx)+log2(1+ nr_bound_MRC/dr_bound_MRC);
            %% ZF SINR
            nr_ZF = pu; invGG=inv(G'*G);
            dr_ZF = invGG(k,k);
            %% RATE FOR ZF
            rate_ZF(mx) = rate_ZF(mx)+log2(1+nr_ZF/dr_ZF);
            %% BOUNDS FOR ZF
            nr_ZF_bound = pu*beta(k)*(M-K);
            bound_ZF(mx) = bound_ZF(mx)+log2(1+ nr_ZF_bound);
        end
    end
end
%%
rate_MRC = rate_MRC/ITER;
bound_MRC = bound_MRC/ITER;
rate_ZF = rate_ZF/ITER;
bound_ZF = bound_ZF/ITER;

figure;
plot(Mv,rate_MRC, 'mx-','LineWidth',2)
hold on
plot(Mv, bound_MRC, 'square', 'LineWidth',2)
hold on
plot(Mv,rate_ZF, 'rx-','LineWidth' ,2)
hold on

plot(Mv, bound_ZF, 'o', 'LineWidth', 2)
grid on
title('Sum-Rate of Massive MIMO with Perfect CSI')
legend('MRC', 'Bound MRC',' ZF','Bound_ZF', 'Location', 'SouthEast');
xlabel('Number of BS Antennas');
ylabel('Uplink Sum Rate (bits/s/Hz)');