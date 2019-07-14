clear all
close all

delta1sq=1;%Expected channel gain of the weaker user;
delta2sq=5;%Expected channel gain of the stronger user;

SNRdB=20:30; SNR=10.^(SNRdB./10);

ITER=100000; %Number of iterations
tildeR_1=1; tildeR_2=1;
R1=2^(tildeR_1)-1;
R2=2^(tildeR_2)-1;

%%
for ix=1:length(SNR)
    ix
    cap_opt=0;
    cap_subopt=0;
    rhos=SNR(ix);
    for k=1:ITER
        h1=sqrt(delta1sq/2).*(randn(1)+1j*randn(1)); %Channel fading coefficient from base station to user 1
        h2=sqrt(delta2sq/2).*(randn(1)+1j*randn(1)); %Channel fading coefficient from base station to user 2
        
        if(abs(h1)>abs(h2))
            h_temp=h1;
            h1=h2;
            h2=h_temp; %Swapping the two values
        end
        
        %To satisfy QoS constraint at two users
        a2_max=(rhos*(abs(h1)^2)-R1)/(rhos*(abs(h1)^2)*(1+R1));
        a2_min=R2/(rhos*(abs(h2)^2));
        
        if(a2_max>=a2_min)
            a2_opt=a2_max;
            a1_opt=1-a2_opt;
        else
            while(a2_max<a2_min)
                h1=sqrt(delta1sq/2).*(randn(1)+1j*randn(1)); %Channel fading coefficient from base station to user 1
                h2=sqrt(delta2sq/2).*(randn(1)+1j*randn(1)); %Channel fading coefficient from base station to user 2
                
                if(abs(h1)>abs(h2))
                    h_temp=h1;
                    h1=h2;
                    h2=h_temp; %Swapping the two values
                end
                a2_max=(rhos*(abs(h1)^2)-R1)/(rhos*(abs(h1)^2)*(1+R1));
                a2_min=R2/(rhos*(abs(h2)^2));
            end
            a2_opt=a2_max;
            a1_opt=1-a2_opt;
        end
        
        %Achievable Sum-Rate With Optimum Power Allocation
        gamma1_u1_subopt=(a1_opt*rhos*abs(h1)^2)/(a2_opt*rhos*abs(h1)^2+1);
        gamma2_u2_opt=a2_opt*rhos*abs(h2)^2;  %For decoding x_2
        cap_opt=cap_opt+log2(1+gamma1_u1_subopt)+log2(1+gamma2_u2_opt);
        
        %Achievable Sum-Rate With Random Power Allocation
        a2=a2_min+rand(1)*(a2_max-a2_min);
        a1=1-a2;
        gamma1_u1_subopt=(a1*rhos*abs(h1)^2)/(a2*rhos*abs(h1)^2+1);
        gamma2_u2_subopt=a2*rhos*abs(h2)^2;
        cap_subopt=cap_subopt+log2(1+gamma1_u1_subopt)+log2(1+gamma2_u2_subopt);
    end
    cap_subopt_avg(ix)=cap_subopt/ITER;
    cap_opt_avg(ix)=cap_opt/ITER;
end
%% Plotting commands
plot(SNRdB,cap_subopt_avg, 'r - s','Linewidth',2.0)
hold on
plot(SNRdB,cap_opt_avg,'k-.o','Linewidth', 2.0)
grid on
legend('Capacity with Random Power Allocation', 'Capacity with Optimum Power Allocation')
xlabel('SNR (dB)')
ylabel('Capacity');
title('Capacity vs SNR')
