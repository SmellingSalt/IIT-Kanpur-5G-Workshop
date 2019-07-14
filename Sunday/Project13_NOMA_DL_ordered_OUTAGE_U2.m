clc
clear all
close all
%Initializations
delta1sq=1;%Expected channel gain of the weaker user;
delta2sq=5;%Expected channel gain of the stronger user;

SNRdB=0:20; SNR=10.^(SNRdB./10);
a1=0.5; a2=0.5;

ITER=100000; %Number of iterations
tildeR_1=0.1; tildeR_2=0.1;

%% Outage Probability Simulation
for ix=1:length(SNR)
    ix
    outage1=0; outage2=0;
    rhos=SNR(ix);
    for k=1:ITER        
        h1=sqrt(delta1sq/2).*(randn(1)+1j*randn(1)); %Channel fading coefficient from base station to user 1        
        h2=sqrt(delta2sq/2).*(randn(1)+1j*randn(1)); %Channel fading coefficient from base station to user 2        
        if(abs(h1)>abs(h2))
            h_temp=h1;
            h1=h2;
            h2=h_temp; %Swapping the two values
        end
        gamma1_u1=(a1*rhos*abs(h1)^2)/(a2*rhos* abs(h1)^2+1);
        gamma1_u2=(a1*rhos*abs(h2)^2)/(a2*rhos*abs(h2)^2+1);
        gamma2_u2=a2*rhos*abs(h2)^2;
        
        %Outage at user 1        
        if(log2(1+gamma1_u1)<tildeR_1)
            outage1=outage1+1;            
        end
        %Outage at user 2
        if((log2(1+gamma1_u2)<tildeR_1) || (log2(1+gamma2_u2)<tildeR_2))            
            outage2=outage2+1;
        end
    end
    Pout1(ix)=outage1/ITER;
    Pout2(ix)=outage2/ITER;
    
    %% Analytical
    R1=2^(tildeR_1)-1;
    R2=2^(tildeR_2)-1;
    phi=max((R2/a2),(R1/(a1-a2*R1)));
    delta3sq=(delta1sq* delta2sq)/(delta1sq+delta2sq);
    Pout1_theory(ix)=1-exp(-R1/(delta3sq*rhos*(a1-a2*R1)));
    Pout2_theory(ix)=1-exp(-phi/(rhos*delta1sq))-exp(-phi/(rhos*delta2sq))+exp(-phi/(rhos*delta3sq));    
end
%%
semilogy(SNRdB, Pout1, 'r', 'Linewidth',2.0)
hold on
semilogy(SNRdB, Pout2, 'k', 'Linewidth',2.0)
hold on
semilogy(SNRdB,Pout1_theory,'s', 'Linewidth',2.0)
hold on
semilogy(SNRdB,Pout2_theory,'o', 'Linewidth',2.0)
grid on
legend('Outage User 1 (Sim.)','Outage User 2 (Sim.)','Outage User 1 (Theory)', 'Outage User 2(Theory)')
xlabel('SNR (dB)')
ylabel('Prob of Outage')
title('Pout vs SNR for Ordered DL NOMA')