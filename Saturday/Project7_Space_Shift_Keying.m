close all
clearvars;
SNRdB=0:2:10;
ITER = 10000;
Nt = 4;
Nr = 4;
bpcu=log2(Nt);
BERopt=zeros(1, length(SNRdB));
BERsubopt = zeros(1, length(SNRdB));
for K=1:length(SNRdB)
    rho = 10^(SNRdB(K)/10);
    for ite = 1:ITER
        disp(ite)
        TxBits=randi([0,1], 1, bpcu);
        H= (1/sqrt(2))*(randn(Nr,Nt)+ 1j*randn(Nr, Nt));
        % H = H* diag(1./(sqrt(sum(abs(H).^2,1))});
        RxNoise =(1/sqrt(2))*(randn(Nr,1)+1j*randn(Nr,1));
        antIndex = 1+bin2dec(num2str(TxBits));
        RxVec = sqrt(rho)*H(:,antIndex) + RxNoise;
        
        % Optimal ML Detector
        ColNorm=sum(abs(H).^2,1);
        MLobj=sqrt(rho)*ColNorm - 2*real(RxVec'*H);
        [minval, minIdx] = min(MLobj);
        DecBits =dec2bin(minIdx-1,log2(Nt))~='0';
        BERopt(1, K)=BERopt(1, K) + sum(DecBits~=TxBits);
        
        % Suboptimal detection
        [maxVal, maxIdx] = max(abs(H'*RxVec));
        DecBits = dec2bin(maxIdx-1,log2(Nt))~='0';
        BERsubopt(1, K)=BERsubopt(1, K) + sum(DecBits~=TxBits);
        
    end
end
%%
BERopt = BERopt/(bpcu*ITER);
BERsubopt = BERsubopt/(bpcu*ITER);
semilogy(SNRdB,BERopt,'--ok', 'LineWidth',2)
hold on
semilogy(SNRdB,BERsubopt, '-sr','LineWidth', 2)
grid
xlabel('SNR (dB)')
ylabel('P_e')
legend('BER Optimal', 'BER Suboptimal')
title('BER performance of Space Shift keying vs SNRI')