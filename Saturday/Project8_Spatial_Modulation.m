clc
close all;
clearvars;
SNRdB = 0:2:10;
ITER = 10000;
Nt=8;
Nr = 8;
M=2;% BPSK
bpcu=log2(Nt*M);
BERopt = zeros(1, length(SNRdB));
BERsubopt=zeros(1, length(SNRdB));
BERtheory = zeros(1,length(SNRdB));
for K=1:length(SNRdB)
    rho =10^(SNRdB(K)/10);
    for ite = 1:ITER
        disp(ite)
        TxBits = randi([0,1], 1, bpcu);
        H= (1/sqrt(2))*(randn(Nr,Nt)+ 1j*randn(Nr,Nt));
        % Normalizing each column for constrained channel condition
        %H = H*diag(1./(sqrt(sum(abs(H).^2,1))));
        RxNoise = (1/sqrt(2))*(randn(Nr,1)+1j*randn(Nr,1));
        iBit = (2*TxBits(1,bpcu)-1);
        antindex = 1+bin2dec(num2str(TxBits(1,1:log2(Nt))));
        RxVec = sqrt(rho)* H(:,antindex)*iBit+RxNoise;
        
        % Optimal ML Detector
        ColNorm = [sum(abs(H).^2,1), sum(abs(H).^2,1)];
        MLobj=sqrt(rho) *ColNorm-2*real(RxVec'*[H,-H]);
        %Least distance is considered
        [minVal, minIdx] = min(MLobj);
        if minIdx <= M*Nt/2
            BitDec = 1;
            antindDec=dec2bin(minIdx-1,log2(Nt))~='0';
        else
            BitDec = 0;
            antindDec = dec2bin(minIdx-Nt-1,log2(Nt))~='0';
        end
        DecBits=[antindDec, BitDec];
        BERopt(1, K) = BERopt(1, K) + sum(DecBits~=TxBits);
        % Suboptimal detection Antenna first
        [maxVal, maxIdx] = max(abs(H'*RxVec));
        ColH = H(:,maxIdx);
        % Decoding of info bits
        bitDec = (real(ColH'*RxVec)>=0);
        antindDec = dec2bin(maxIdx-1,log2(Nt))~='0';
        DecBits = [antindDec, bitDec];
        
        BERsubopt(1, K) = BERsubopt(1, K) + sum(DecBits ~=TxBits);
    end
    
    BERtheory(1, K) = BERtheorySM(SNRdB(K), Nt, Nr, M);
end
%%
BERopt = BERopt/(bpcu*ITER);
BERsubopt = BERsubopt/(bpcu*ITER);
semilogy(SNRdB,BERopt,'--ok','Linewidth',2)
hold on
semilogy(SNRdB,BERsubopt,'-*r', 'LineWidth',2)
semilogy(SNRdB,BERtheory,'-.^g','LineWidth',2)
xlabel('SNR (dB)')
ylabel('P_e')
legend('BER Optimal', 'BER Sub-optimal', 'BER Theory')
title('BER performance of Spatial Modulation vs SNR')
grid on;
axis tight;