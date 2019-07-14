close all; clear all; rng('shuffle'); clc
% simulation parameters
t=32; r=32; %% Number of Tx/Rx Antennas
numRF=6; %% Number of RF Chains
G=64; %% Grid Size
L=8; %% Sparsity level
Ns=6; %% Number of data streams
ITER=100; %% Number of iterations
% Initializations
H = zeros(r,t); %% channel matrix
SNRdB=-5:1:5;
C_HYB=zeros(length(SNRdB),1); %% Capacity of Hybrid MIMO
C_MIMO= zeros(length(SNRdB),1);%% Capacity of Conventional MIMO
% G-quantized Array response matrix
A_T =zeros(t,G);

for I=1:G
    dirCos=2/G*(I-1)-1;
    for K=1:t
        A_T(K,I)= 1/sqrt(t)*exp(-1j*pi*(K-1)*dirCos);
    end
end
A_R=A_T; %% for simplicity
for iter1=1:ITER
    fprintf("First loop, iteration %d\n",iter1);
    % generation of channel
    A_T_genie=zeros(t,L); A_R_genie = zeros(r,L);
    % generate AoD/AoAs uniformly in grid
    AoDlist = randperm(G,L); AoAlist=randperm(G,L);
    % tx/rx array response matrix corr. to true AOD/AOAS
    A_T_genie=A_T(:,AoDlist); A_R_genie=A_R(:,AoAlist);
    % generate channel gain
    chGain=1/sqrt(2)*(randn(L,1)+1j*randn(L,1));
    % generate channel matrix
    H=sqrt(t*r/L)*A_R_genie*diag(chGain)*A_T_genie';
    % SVD of H
    [U,S,V]=svd(H);
    % optimal unconstrained precoder
    Fopt=V(:,1:Ns);
    % OMP based transmit hybrid precoder design
    [FBB, FRF] = SOMP_mmW_precoder(Fopt, A_T, eye(t), numRF);
    % NORM to meet power constraint
    FBB_NORM=sqrt(Ns)/(norm(FRF*FBB, 'fro'))*FBB; %% normalized precoder
end
for i_snr=1:length(SNRdB)
    %% MMSE COMBINER
    np = 10^(-SNRdB(i_snr)/10); %% noise power for signal power 1
    % optimal unconstrained MMSE combiner for unconstrained precoder
    % Standard mmse combiner equation
    Wmmse_opt = H*Fopt*inv(Fopt'*H'*H* Fopt + np*Ns*eye(Ns));
    % Capacity of optimal unconstrained MIMO precoder/combiner
    %Standard formula for MIMO capacity
    C_MIMO(i_snr) = C_MIMO(i_snr) + mimo_capacity(Wmmse_opt'*H*Fopt, 1/Ns*eye(Ns),np*Wmmse_opt'*Wmmse_opt);
    % optimal unconstrained MMSE combiner for hybrid precoder
    %Covariance of the recieved vector
    Ryy=1/Ns*H*FRF*FBB_NORM*FBB_NORM'*FRF'*H' +np*eye(r);
    
    %% DESIGNING HYBRID COMBINER
    Wmmse_Hyb=H*FRF*FBB_NORM*inv((FBB_NORM')*(FRF')*(H')*H*FRF*FBB_NORM+np*Ns*eye(Ns));
    % OMP based receive hybrid combiner design
    [WBB, WRF] = SOMP_mmW_precoder(Wmmse_Hyb, A_R, Ryy, numRF);
    % Capacity of OMP-based precoder/combiner
    C_HYB(i_snr)=C_HYB(i_snr) + mimo_capacity((WBB')* (WRF')*H*FRF*FBB_NORM, 1/Ns*eye(Ns),np*(WBB')*(WRF')*WRF*WBB);
end
C_MIMO=C_MIMO/ITER; C_HYB=C_HYB/ITER;
plot(SNRdB,abs(C_MIMO),'b','linewidth',3.0);
hold on;
plot(SNRdB,abs(C_HYB),'m -.','linewidth',3.0);
grid on; axis tight;
xlabel('SNR(dB)'); ylabel('Capacity (b/s/Hz)');
legend('Conventional MIMO','Hybrid MIMO'); title('Capacity vs SNR');
%end of the code