clc; close all; clear all; rng('shuffle');
% simulation parameters
t=32; r=32; %% Number of Tx/Rx Antennas
numRF = 8; %% Number of RF Chains
N_Beam=24; %% Number of Pilot Symbols
G=32; %% Grid Size
ITER=10; %% Number of iterations
L=5; %% Sparsity level
% initializations
omp_thrld=1; %Threshold for omp stopping
kp=zeros(t*r,L); %This is used to create the effective channel (using kroneker product)
SNRdB=10:10:50;
mseOMP=zeros(length(SNRdB),1); %Finding the mean square error using OMP
mseGenie = zeros(length(SNRdB),1);%Finding the actual mean square error using the defined H and Q
%% CONSTRUCTING THE DICTIONARY MATRIX
% G-quantized Array response matrix
A_T=zeros(t,G); A_R = zeros(r,G);
for I=1:G
    dirCos =2/G*(I-1)-1;
    for K=1:t
        A_T(K,I)=1/sqrt(t)*exp(-1j*pi*(K-1)*dirCos);
    end
end
%using the same dictionary at both tx and rx;
A_R=A_T; %% For simplicity
%% PRECODER AND COMBINER CREATION
% load the tx/rx precoder/combiner matrices
path=pwd;
pwd=[pwd '\Matlab code for download\mmWave_matrices.mat'];
load(pwd);
%% EFFECTIVE CHANNEL
Q = kron((FBB.')*(FRF.'),(WBB')*(WRF'));
%% SIMULATION
for ix=1:ITER
    disp(ix);
    %% INITIALISNG THE ACTUAL ANGLE OF ARRIVAL AND DEPARTURES
    A_T_genie=[]; A_R_genie=[];
    H=zeros(r,t); %% Channel Matrix
    % generate AoD/AoA uniformly in grid
    %% CHOOSING 'L' PATHS RANDOMLY FROM THE DICTIONARY
    for I=1:L
        ix1=randi([1,G]); %% AoD index
        ix2=randi([1,G]); %% AOA index
        % generate channel gain
        chGain = 1/sqrt(2)*(randn(1,1) + 1j*randn(1,1));
        % CHANNEL MATRIX
        H=H+sqrt(t*r/L)*chGain *A_R(:,ix2)*(A_T(:,ix1))';
        % obtain tx/rx array response matrix corresponding to true Asos/AiDs
        A_T_genie=[A_T_genie, A_T(:,ix1)];
        A_R_genie = [A_R_genie, A_R(:,ix2)];
        %USED TO CREATE EFFECTIVE CHANNEL
        kp(:,I) = kron(conj(A_T(:,ix1)), A_R(:,ix2)); 
    end
    %% NOISE VECTOR VARYING WITH SNR
    % generate the noise Vector
    ChNoise = 1/sqrt(2)*(randn(N_Beam*N_Beam,1)+ 1j*randn(N_Beam*N_Beam,1));
    %% ESTIMATING CHANNEL AND ERRORS USING OMP (NOT KNOWING ORIGINAL CHANNEL) AND LS (KNOWING ORIGINAL CHANNEL)
    for i_SNR=1:length(SNRdB)
        snr = 10^(SNRdB(i_SNR)/10);
        % Measurement Vector
        y=sqrt(snr)*Q*H(:)+ ChNoise(:);
        % equivalent dictionary matrix for CS-problem
        %% EFFECTIVE CHANNEL WITH THE REQUIRED POWER SCALING
        Qbar=sqrt(snr)*Q*(kron(conj(A_T),A_R));
        %% OMP TO ESTIMATE CHANNEL
        %GETTING THE SPARSE AND DIAGONAL hb
        h_b_omp=OMP_mmWave_Est(y,Qbar,omp_thrld); %% estimate of beamspace channel
        %GETTING THE ACTUAL H USING AT,hb and Ar
        % error metric omp
        H_omp = A_R*(reshape(h_b_omp,r,t))*A_T'; %% estimate of channel matrix
        mseOMP(i_SNR) = mseOMP(i_SNR) + ((norm(H-H_omp,'fro'))^2/(t*r));
        
        %% LEAST SQUARE ESTIMATION
        % Computing the pseudo inverse on the known At, Ar and sparse
        % hb(kp)
        Q_ORACLE = sqrt(snr)*Q*kp;
        chGainEst=pinv(Q_ORACLE)*y;
        % error metric for the LS estimated channel
        H_genie=A_R_genie*diag(chGainEst)*A_T_genie';
        mseGenie(i_SNR) = mseGenie(i_SNR) + ((norm(H-H_genie, 'fro'))^2/(t*r));
        
    end
end
%%
mseOMP= mseOMP/ITER; mseGenie = mseGenie/ITER;
% plots
semilogy(SNRdB,mseOMP,'g *-','linewidth',3.0);
hold on;
semilogy(SNRdB, mseGenie,'m o-','linewidth',3.0);
axis tight; grid on;
xlabel('SNRdB'); ylabel('Normalized MSE');
legend('OMP','ORACLE LS'); title('MSE vs SNRdB');
%end of the code