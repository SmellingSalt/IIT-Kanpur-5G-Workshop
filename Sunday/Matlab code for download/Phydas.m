function p_k = Phydas(L_f,N)
H1=0.971960;H2=sqrt(2)/2;H3=0.235147;
fh=1+2*(H1+H2+H3);
hef(1:L_f+1)=0;
for i=0:L_f
   hef(1+i)=1-2*H1*cos(pi*i/(2*N))+2*H2*cos(pi*i/N)-2*H3*cos(pi*i*3/(2*N));
end
hef=hef/fh;
p_k=hef/norm(hef);