%Implementation of system model of figure 1 of transaction analysis and
% design of OFDM/OQAM Systems"
% Follow the block diagram to understand the code

clear all; close all;clc;
ITER = 100;
nSym=50;
N =64; M = N/2; % Number of subcarriers and upsampling factor
L_f= 4*N; % Filter order in the multiple of N
L_h=6; % number of channel taps
L_s=L_f+nSym*M+L_h-1; %length of the received signal
Beta=0; Alpha = L_f/M; % defined according to equation 26 of paper L f = Alpha*M-Beta
rf = 1; % Roll off factor of RRC filter
EsN0dB = 0:3:21;
EsN0=10.^(EsN0dB/10);
BER = zeros(1,length(EsN0dB));

% Prototype filter
%p_k=rcosdesign(rf,L_f/N,N);%generating the impulse response of raised cosine prototype filter
p_k= Phydas(L_f,N);%generating the impulse response of PHYDYAS filter

% Transmit and receive Filterbanks
Tr_FB=[]; Rx_FB=[];
T=(L_f-M)/2;
R=(L_f+M)/2;
%% GENERATING MODULATED FILTER BANK RESPONSE
for m=0:N-1
    Tr_FB = vertcat(Tr_FB,exp(-11*2*pi*(1/N)*m*T)*p_k.*exp(11*2*pi*m*(1/N)*[0:L_f]));
    Rx_FB = vertcat(Rx_FB,exp(-1i*2*pi*(1/N)*m*R)*p_k.*exp(1i*2*pi*m*(1/N)*[0:L_f]));
end


for jj=1:ITER
    jj

    h = 1/sqrt(2)*sqrt(1/L_h)*(randn(1,L_h) + 1i*randn(1,L_h));
    %cHANEL FREQ RESP
    cfr = fft(h,N).'; % Channel frequency response
    nt_Base=randn(1,L_s) + 1i*randn(1,L_s);
    BitFrame=2*randi([0,1],N,nSym)-1;
    nErr=0;
    %BPSK Symbols
    for kk=1:length(EsN0dB)
        %% Transmitter
        Phi_mn = repmat(exp(1i*(pi/2)* [0:nSym-1]),N,1);%calculating exp(j*n*pi*/2)
        SymMatrix=BitFrame.*Phi_mn;
        UpSym=upsample(SymMatrix.',M).';
        TxOut=[];
        for i=1:N
            TxOut=vertcat(TxOut,conv(UpSym(1,:),Tr_FB(1,:)));
        end
        TxSignal=sum(TxOut);
        %% receiver
        RxSignal=conv(h, TxSignal) + sqrt(1/(EsN0(kk)))*nt_Base; % generating Noise
        RxOut = [];
        for xx=1:N
            RxOut=vertcat(RxOut,conv(RxSignal,Rx_FB(xx,:)));
        end
        
        DownSym=downsample(RxOut.',M).';
        RX_Phi_mn = repmat(exp(-(i)*(pi/2)*([0:length(DownSym(1,:))-1]-Alpha)),N,1);
        DownSym = DownSym.*RX_Phi_mn;
        % Constant channel over the coherence bandwidth for all subcarriers
        cfr_mat = repmat(cfr,1,nSym);
        %Equalising the CFR
        EqSym=DownSym(:, Alpha+1:Alpha+nSym)./cfr_mat;
        %Only real part is orthogonal
        RXOQAM=real(EqSym);
        %Quantizing the information
        iHat=2*(RXOQAM>=0)-1;
        %Counting errors
        nErr=sum(sum(BitFrame-iHat~=0));
        %BER for each SNR
        BER(kk)= BER(kk)+ nErr/(N*nSym*ITER);
    end
end
%% plotting
Theory_BER = 0.5*(1-sqrt(EsN0./(EsN0+2)));
figure
semilogy(EsN0dB,Theory_BER,'mx-','LineWidth',2);
hold on
semilogy(EsN0dB,BER, 'bs-', 'LineWidth',2);
grid on
legend('theory', 'simulation');
xlabel('SNR, dB')
ylabel('BER')
axis tight;