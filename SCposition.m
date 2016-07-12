function [r,v] = SCposition(acend_node,i,w,a,e,t)

mue=3.986*10^5;

%Convert to radians
acend_node=acend_node*pi/180;
i=i*pi/180;
w=w*pi/180;
T=2*pi*sqrt(a^3/mue);
M=2*pi*t/T;

E_guess = (M.*(1-sin(M+e))+(M+e).*sin(M))./(1+sin(M)-sin(M+e)); 

functE = @(E)sqrt(a.^3./mue).*(E-e.*sin(E))-t; 

E=fzero(functE,E_guess);

r_pfmag = a.*(1-e.*cos(E));

R= a.*[(cos(E)-e),(sqrt(1-e.^2).*sin(E)),0];
V = sqrt(mue.*a)./r_pfmag.*[-sin(E), sqrt(1-e.^2).*cos(E),0];

%Transformation Matrix
co=cos(acend_node);
so=sin(acend_node);
cw=cos(w);
sw=sin(w);
ci=cos(i);
si=sin(i);

C=[co*cw-so*ci*sw so*cw+co*ci*sw si*sw;...
    -co*sw-so*ci*cw -so*sw+co*ci*cw si*cw;...
    so*si -co*si ci];
C_inv=inv(C);

%Geocentric-Equatorial Coord.
r=C_inv*R';
v=C_inv*V';