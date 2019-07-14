clc; close all; clear all; rng('shuffle');
% simulation parameters
t=32; r=32; %% Number of Tx/Rx Antennas
numBER = 8; %% Number of RF Chains
N_Beam=24; %% Number of Pilot Symbols
G=32; %% Grid Size
ITER=10; %% Number of iterations
L=5; %% Sparsity level
% initializations
omp_thrld=1;
kp zeros(t*r,L);
SNRdB 10:10:50;
mse cOMP zeros(length(SNR dB),1);
massGenie = zeros(length(SNRdB),1);
% G-quantized Array response matrix
A_T=zeros(t,G); A_R = zeros(r,G);
for I=1:G
    dirCos =2/G*(1-1)-1;
    for K=1:t
        A_T(K,I)=1/sqrt(t)*exp(-1j*pi*(K-1)*dirCos);
    end
end
A_R=A_T; %% For simplicity
% load the tx/rx precoder/combiner matrices
load('mmWave matrices');
Q = kron((FBB.')*(FRF.'),(WBB')*(WRF'));
for ix=1:ITER
    disp(ix);
    A_T_genie=[]; A_R_genie=[];
    H=zeros(r,t); %% Channel Matrix
    % generate AoD/AoA uniformly in grid
    for I=1:L
        ix1=randi([1,G]); %% AoD index
        ix2=randi([1,G]); %% AOA index
        % generate channel gain
        chGain = 1/sqrt(2)*(randn(1,1) + 1j*randn(1,1));
        % obtain channel matrix
        H=H+sqrt(t*r/L)*chGain *A_R(ix2)*(A_T(:,ix1))';
        % obtain tx/rx array response matrix corresponding to true Asos/AiDs
        A_T_genie=[A_T_genie, A_T(:,ix1)];
        A_R_genie = [A_R_genie, A_R(:,ix2)];
        kp(:,I) = kron(conj(A_T(:,ix1)), A_R(:,ix2));
    end
    % generate the noise Vector
    ChNoise = 1/sqrt(2)*(randn(N_Beam*N_Beam,1)+ 1j*randn(N_Beam*N_Beam,1));
    for i_SNR=1:length(SNRdB)
        snr = 10^(SNRdB(i_SNR)/10);
        % Measurement Vector
        y=sqrt(snr)*Q*H(:)+ ChNoise(:);
        % equivalent dictionary matrix for CS-problem
        Qbar=sqrt(snr)*Q*(kron(conj(A_T),A_R));
        % call OMP function
        h_b_omp=OMP_mmWave_Est(y,Qbar,omp_thrld); %% estimate of beamspace channel
        % error metric omp
        Homp = A_R*(reshape(h_b_omp,r,t))*A_T'; %% estimate of channel matrix
        mseOMP(i_SNR) = museOMP(i_SNR) + ((norm(H-H_omp,'fro'))^2/(t*r));
        % ORACLE LS
        Q_ORACLE = sqrt(snr)*Q*kp;
        chGainEst=pinv(Q_ORACLE)*y;
        % error metric genie
        H_genie=A_R_genie*diag(chGainEst)*A_T_genie';
        mseGenie(i_SNR) = mseGenie(i_SNR) + ((normal(H-H_genie, 'fro'))^2/(t*r));
        
    end
    
end

museUM = museUM/ITER; massGenie = massGenie/ITER;
% plots
semilogy(SNRdB,mseOMP,'g *-','linewidth',3.0);
hold on;
semilogy(SNRdB, mseGenie,'m o-','linewidth',3.0);
axis tight; grid on;
xlabel('SNRdB'); ylabel('Normalized MSE');
legend('OMP','ORACLE LS'); title('MSE vs SNRdB');
%end of the code