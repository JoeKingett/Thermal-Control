function [qtot,Fsce,qs,qa,qe] = ext_heat_flux_loop(r_sc,r_sun,n,Gs,AF,abso_sol,emis_ir,refl_cb,r_mag)

% This function calculates the external (direct solar, qs, albedo radiation, 
% qa and Earth IR) heat flux (W/nT2) absorbed by each input surface. The 
% output is a 1 x 3 x nsides vector, where nsides is defined below. Each 
% page of the output vector corresponds to one of the sides of the 
% spacecraft. The 3 columns correspond to[qs,qa,qe] respectively. The total 
% heat flux would be sum(output),and the total heat flow absorbed by the 
% surface (W) would be sum (output)*area of the surface.

% Inputs:
% r_sc = surf, position geocentric-equitorial coords [rx,ry,rz]) (Km)
% r_sun = position of the sun in geocen coords [rx,ry,rz] (km)
% n = normal vector of surface in geocen coords [nx,ny,z]
% Gs = solar constant (W/nT2)
% AF = Albedo Factor
% Te = Effective BB temperature of Earth (k)
% abso_sol = solar absorbtivity of surface (unitless)
% emis_ir = Infrared emissivity of surface (unitless)

% Constants
sigma = 5.6704e-8;                   %Steffan Boltzmann const (W/(m~2 K))
Re=6378.14;                          %radius of Earth (km)
sizen = size(n);                %dimensions of the surfaces' normals matrix
nsides = 1;                     %the # of pages of n = # of surfaces
qtot = zeros(1,3,nsides);       %pre-alocate, to be filled in loop

for k=1:nsides
h = norm(r_sc)-Re; %SC altitude
gammar=acos(-sum(n(:,:,k).*r_sc)./(norm(n(:,:,k)).*norm(r_sc))); 
                            %angle between n and r_sc in radians
gamma = gammar.*180./pi; 
                            %angle between n and r_sc in deg
Fsce = view_factor_scearth(h,gamma);
                            %view factor between Earth and surface
S = r_sun-r_sc; 
                            %line of sight vector between sc and sun (km)
psir = acos(sum(n(:,:,k).*S)./(norm(n(:,:,k)).*norm(S)));
                %angle between line of sight vector and surface normal vector (radians)
thetar = acos(sum(r_sun.*r_sc)./(norm(r_sun).*norm(r_sc)));
                  %solar reflection angle off earth (radians)

%%Direct Solar radiation absorbed by sc per sqr meter%%
qsa = Gs.*abso_sol.*cos(psir);
qsr = refl_cb.*qsa;
qs = qsa-qsr;

%%Reflected Solar (albedo) radiation absorbed by sc W/nT2%%
qaa = Gs.*AF.*Fsce.*abso_sol.*cos(thetar);
qar = refl_cb.*qaa;
qa= qaa-qar;
%%Direct Earth IR radiation absorbed by sc W/m~2%%
%qe = 237*(Re/(r_mag))^2*emis_ir*Fsce; %ir abs.=ir emiss.
qe=237*sigma*emis_ir*Fsce; %ir abs.=ir emiss.

if insun(r_sc,r_sun)==0;    %If in ecclipse:
qs = 0;                         %overwrite qs to zero
qa = 0;                         %overwrite qa to zero
end

if psir >= pi/2;            %If surface pointing away from sun:
qs = 0;                             %overwrite qs to zero
end

if qa < 0                   %albedo goes to zero for theta>pi/2
qa = 0;
end


%output environmental radiation abosrbed by surface (W/nT2)
qtot(:,:,k)=[qs,qa,qe];

end
end




