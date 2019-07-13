clear;
N = 20; M = 10;
Phi =2*randi ([0,1],M,N) -1;
f = [0 0 1.6 0 0 0 0 0 5.5 0 0 0 0 0 0 0 0 0 1.2 0]';
y = Phi*f;% // Measurement equation
Phi_hat=[];
res=y;
nz= [];
f_omp = zeros (20,1);
for ix = 1:3
    [cor, col]= max(abs (Phi'*res));
    nz =[nz;col];
    Phi_hat=[Phi_hat, Phi(:,col)];
    hat_f= pinv(Phi_hat)*y;
    nz =[nz;col];
    Phi_hat=[Phi_hat,Phi(:,col)];
    hat_f=pinv(Phi_hat)*y;
    res=y - Phi_hat*hat_f;
end
f_omp(nz)=hat_f;
f_min=(Phi')*inv(Phi*Phi')*y;
scatter ([1:N], f,'b o', 'linewidth', 1.0);
% W Scatter plot of original signal
hold on;
grid;
scatter([1:N],f_omp, 'm s','linewidth',1.0);
figure;
scatter ([1:N],f_min, 'k s','linewidth',1.0);