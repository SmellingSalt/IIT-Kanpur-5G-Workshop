function p_k = rrc(rf,N,T_s)
T_int = T_s/N; %sampling interval
Tsi=-2*T_s:T_int:2*T_s; % sampling instants
p_k=zeros(1,length(Tsi));
    for i=1:1:length(Tsi)
        if Tsi(i)==0
            p_k(i)= (1-rf)+4*rf/pi;
        else if Tsi(i)==1/(4*rf) || Tsi(i)==-1/(4*rf)
               p_k(i)=rf/sqrt(2)*((1+2/pi)*sin(pi/(4*rf))+(1-2/pi)*cos(pi/(4*rf)));
              else
                p_k(i) = (sin(pi*Tsi(i)*(1-rf))+4*rf*Tsi(i).*cos(pi*Tsi(i)*(1+rf)))./(pi*Tsi(i).*(1-(4*rf*Tsi(i)).^2));
             end
        end
    end
p_k = p_k/norm(p_k);
