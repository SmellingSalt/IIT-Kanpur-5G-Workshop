%% GRAPH NOT MATCHING %%
clc;
ITER = 1000;
K=10;
Mv=20:30:500;
EudB=10;
Eu=10^(EudB/10);
nTsym=K; %Number of pilots
rateMRC=zeros(1,length(Mv));
boundMRC=zeros(1,length(Mv));
rateZF=zeros(1,length(Mv));
boundZF=zeros(1,length(Mv));

for ite=1:ITER
    disp(ite)
    D=Dmatrix(K); ubeta = diag(D);
    for mx=1:length(Mv)
        
        M=Mv(mx);
        %pu = Eu;
        %pu = Eu/M;
        pu = Eu/sqrt(M);
        Pp=nTsym*pu; %Pilot power 
        Phi=sqrt(1/nTsym)*dftmtx(nTsym); %Pilot matrix
        pilotMtx =sqrt(Pp)*Phi;
        %% CHANNEL
        H = sqrt(1/2)*(randn(Mv(mx),K)+1i*randn(Mv(mx),K));
        G=H*sqrt(D);
        %% 
        N = sqrt(1/2)*(rand(M,nTsym)+1i*randn(M,nTsym));
        RxBlk = G*(pilotMtx.')+N; %recieved pilot matrix
        Ghat =sqrt(1/Pp)*RxBlk*conj(Phi)*inv(1/Pp*inv(D)+eye(K)); %Estimated pilot
        invGG =inv(Ghat'* Ghat);
        for k=1:K
            gkHat=Ghat(:,k);%Extract kth column of the estimated channel
            MUImrc=0; boundUIMmrc = 0;
            errZF = 0; errMRC=0;
            DrMRC=norm(gkHat)^2;
            %Calculate rate of each user and sum it up
            for iu=1:K
                DrErr = (nTsym*pu*D(iu,iu))+1;
                DrMRC=DrMRC +(ubeta(iu)/DrErr)*pu* norm(gkHat)^2;
                errZF=errZF+pu*D(iu,iu)/DrErr;
                if (iu~=k)
                    DrMRC=DrMRC+pu*abs(gkHat'*Ghat(:,iu))^2;
                end
            end
            NrMRC=pu*norm(gkHat)^4;
            rateMRC(mx)=rateMRC(mx)+ log2(1+NrMRC/DrMRC);
            
            NrBoundMRC = nTsym*pu*(M-1)*ubeta(k)^2;
            DrBoundMRC = (nTsym*pu*ubeta(k)+1)* (sum(ubeta)-ubeta(k))+(nTsym+1)* ubeta(k)+(1/pu);
            boundMRC(mx) = boundMRC(mx)+ log2(1+NrBoundMRC/DrBoundMRC);
            
            NrZF=pu;
            DrZF=(errZF+1)*invGG(k,k);
            rateZF(mx) = rateZF(mx)+log2(1+NrZF/DrZF);
            
            NrBoundZF=nTsym*pu*pu*(M-K)*ubeta(k)^2;
            DrBoundZF = (nTsym*pu*ubeta(k)+1)*errZF+(nTsym*pu*ubeta(k))+1;
            boundZF(mx)=boundZF(mx)+ log2(1+NrBoundZF/DrBoundZF);
        end
    end
end
%% PLOT
rateMRC = rateMRC/ITER;
boundMRC = boundMRC/ITER;
rateZF = rateZF/ITER;
boundZF=boundZF/ITER;
figure;
plot(Mv,rateMRC,'mx-', 'LineWidth',2)
hold on
plot(Mv,boundMRC, 'square', 'LineWidth', 2)
hold on
plot(Mv,rateZF, 'rx-', 'LineWidth', 2)
hold on
plot(Mv, boundZF, 'o', 'LineWidth', 2)
grid on
title('Sum-Rate of Massive MIMO with Imperfect CSI')
legend('MRC', 'Bound MRC', 'ZF', 'Bound ZF', 'Location', 'SouthEast');
xlabel('Number of BS Antennas')
ylabel('Uplink Sum Rate (bits/s/Hz)')