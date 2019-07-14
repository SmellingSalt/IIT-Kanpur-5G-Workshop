clear all
close all

delta1sq=5; %Expected channel gain of the stronger user
delta2sq=1;%Expected channel gain of the weaker user

SNRdB=0:2:30; SNR=10.^(SNRdB./10);
a1=0.9; a2=0.1;
ITER=500000; %Number of iterations
tildeR_1=1; tildeR_2=1;
%% Outage Probability Simulation
for ix=1:length(SNR)
    fprintf("Computing error for SNR %d dB \n",ix-1);
    outage1=0; outage2=0;
    rhos=SNR(ix);
    for k=1:ITER        
        h1=sqrt(delta1sq/2).*(randn(1)+1j*randn(1)); %Channel fading coefficient from base station to user 1        
        h2=sqrt(delta2sq/2).*(randn(1)+1j*randn(1)); %Channel fading coefficient from base station to user 2
       
        gamma1=(a1*rhos*abs(h1)^2)/(a2*rhos*abs(h2)^2+1); %For decoding x 1 at BS
        gamma2=a2*rhos*abs(h2)^2;      %For decoding x_2 at BS
        %Outage at user 1
        if(log2(1+gamma1)<tildeR_1)           
            outage1=outage1+1;            
        end       
        %Outage at user 2
        if((log2(1+gamma1)<tildeR_1) || (log2(1+gamma2)<tildeR_2))            
            outage2=outage2+1;
        end
    end
    Pout1(ix)=outage1/ITER;
    Pout2(ix)=outage2/ITER;
    
    %% Analytical
    R1=2^(tildeR_1)-1;
    R2=2^(tildeR_2)-1;
    phi2=(R1/(a1*delta1sq))+(R2/(a2*delta2sq))+((R1*R2)/(a1*delta1sq));
    phi3=1+((R1*a2*delta2sq)/(a1* delta1sq));
    Pout1_theory(ix)=1-((exp(-R1/(a1*rhos*delta1sq)))/(1+((R1*a2*delta2sq)/(a1*delta1sq))));
    Pout2_theory(ix)=1-((exp(-phi2/rhos)/(phi3)));    
end

%%

semilogy(SNRdB, Pout1,'r','Linewidth',2)
hold on
semilogy(SNRdB, Pout2,'k', 'Linewidth',2)
hold on
semilogy(SNRdB, Pout1_theory,'s', 'Linewidth',2);
hold on
semilogy(SNRdB, Pout2_theory,'o', 'Linewidth',2)

grid on
legend('Outage User 1 (Sim.)', 'Outage User 2 (Sim.)','Outage User 1 (Theory)','Outage User 2(Theory)')
xlabel('SNR (dB)')
ylabel('Prob of Outage')
title('Pout vs SNR for UL fixed NOMA')

