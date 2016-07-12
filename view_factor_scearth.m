function Fe = view_factor_scearth(h,gamma)
%This function calculates the view factor (also known as configuration
%factor, geometry factor) Fdl-2, from an infinitesimal surface (dl)
%to the Earth (2). Output, Fe, is unitless.

%Inputs:
%h = a l t i t u d e of the surface above E a r t h ' s surface (km)
%gamma = angle between normal vector and nadir vector (degrees)

gammar = gamma.*pi./180; %gamma in radians
Re=6378.14; %Radius of the Earth (km)
rsc=Re+h; %distance of spacecraft from center of Earth (km)
H = rsc./Re;
phi_m = asin(1/H);
b = sqrt(H.^2-1);
%if full Earth is visible to the plate
if gammar <= pi./2-phi_m;
Fe = cos(gammar)./H.^2;
%if part of the Earth is visisible to the plate
elseif gammar > pi./2-phi_m && gammar <= pi./2+phi_m;
tl = 1./2.*asin(b./(H.*sin(gammar)));
t2 = 1./(2.*H.^2).*(cos(gammar).*acos(-b.*cot(gammar))...
-b.*sqrt(1-H.^2.*(cos(gammar))^2));
Fe = 2./pi.*(pi./4-tl+t2);
%if none of the Earth is visible to the plate
else;
Fe=0;
end
end






