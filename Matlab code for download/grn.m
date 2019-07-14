clear all;
r = 6;  %Indices for conversion matrix
N = 64; %Number of Subcarrier
L = 4;  %Oversampling factor
len = L*N ;
Grn = [] ;
for n = 0:2  %Indices for conversion matrix
t = [1 1 1 1 1 1 1 1 ; 1 1 -1 -1 1i 1i -1i -1i ;1 -1 1 -1 1 -1 1 -1 ; -1 1 1 -1 -1i 1i 1i -1i] ;
a = len/(2^(n+2)) ;

g = [t(1,r) zeros(1,a-1),t(2,r),zeros(1,a-1),t(3,r),zeros(1,a-1),t(4,r),zeros(1,a-1),zeros(1,len - (len/(2^n)))].' ;
G = [] ;
for i = 1:len
  G = [G , circshift(g,(i-1))] ;  
end
 Grn = [Grn ; G] ;
end
Grn = 0.5*Grn ;
