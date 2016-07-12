function insun = insun(r_sc,r_sun)

%This function determines whether a spacecraft is in sunlight
%(returning Dor in eclipse (returning 0) given the inputs:

%r_sc = spacecraft position in geocentric coords(km)
%r_sun = sun position in geocentric coords (km)

Re = 6378; %radius of Earth (km)

theta1 = acos(Re./norm(r_sc)); %angle (see note book p29a) (rads)
theta2 = acos(Re./norm(r_sun)); %angle (see note book p29a) (rads)

psi = acos(sum(r_sc.*r_sun)./(norm(r_sc).*norm(r_sun)));
    %angle between sc position vector and sun position vector (rad)\

%if psi is <= thetal+theta2, its in sun, otherwise its in eclipse
if psi >= theta1+theta2;
    insun=0;
else
    insun=1;
end
end