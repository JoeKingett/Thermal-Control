%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-----------------------------------------------------------------------%%
%%----------------------------THERMAL MODEL------------------------------%%
%%-----------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear,clc

%----Orbital Paramters----%
mue=3.986*10^5;                 %gravitational paramater
Re=6348.14;                     %radius of Earth(km)
acend_node=1.422;                   %Right accension of accending node(deg)
i=51.64;                            %inclination(deg)
w=28.6432;                      %Argument of perigee(deg)
a=6782.612778;                  %semi-major axis(km)
e=.0000798;                     %eccentricity(deg)
T=2*pi*sqrt(a^3/mue);           %Orbital Period(s)
delta_t=10;                     %time step (s)
t_final=T*5;                    %final time (s)

%----Constants----%
Gs=1371;                        %Solar Constant(W/m2)
AF=0.367;                       %albedo Factor 
sigma = 5.6704e-8;              %Stefan-Boltzmann constant (W/m2K4)

%----CubeSat Properties----%
ntot=6;                                                         %number of nodes
nsides=6;                                                       %number of sides (outside panels)
A=[0.01 0.01 0.01 0.01 0.01 0.01];                              %surface area per node(m2)
th=[0.006 0.006 0.006 0.006 0.006 0.006];                       %thickness per node (m)
V=A.*th;                                                        %node volumes (m^3)
%Properties of Aluminum
rho=[2712 2712 2712 2712 2712 2712];                            %density (kg/m3)
cp=[961 961 961 961 961 961 ];                                  %specific heat (J/kg/K)
C=rho.*V.*cp;                                                   %Node Thermal Capacities
abso_sol=[ 0.70 0.70 0.70 0.70 0.70 0.70];                      %Solar Absorbance 
emis_ir=[0.79 0.79 0.79 0.79 0.79 0.79];                        %Infrared Emittance
refl_cb = 1-abso_sol;                                           %relfectivity of cube surface

%----Inital Values----%
temp_init=[ 300 300 300 300 300 300];                       %Inital temperature (K)
rad_vf= [   0    0.45   0.2    0.2    0.2    0.2;
          0.45    0      0.2    0.2    0.2    0.2;
           0.2     0.2     0   0.45   0.2    0.2;
           0.2     0.2   0.45    0     0.2    0.2;
           0.2     0.2     0.2    0.2     0    0.45;
           0.2     0.2     0.2    0.2   0.45    0  ];     %Radiation view factor between nodes
       
conductance=[ 0 0 0 0 0 0;                                  %Conductance between nodes (W/K)
              0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0];
          
inthl= [1 1 1 1 1 1];                                       %internal heat load

       
%----Sun position----%
au=149597900;                                              %conversion factor (km/au)
sol_pos=[0.984 0.9888 0.9962 1.005 1.0122 1.0163 1.0161 1.0116 1.0039 0.9954 0.9878 0.937];
month=6;                                                   %month of the year
rho_sun= sol_pos(month)*au;
phi_sun=[113.5 105.67 97.83 90 82.17 74.33 66.5 74.33 82.17 90 97.83 105.67];
theta_sun=[180 210 240 270 300 330 0 30 60 90 120 150];
r_sun=[rho_sun*cosd(phi_sun(month))*cosd(theta_sun(month)) rho_sun*cosd(phi_sun(month))*sind(theta_sun(month))  rho_sun*sind(phi_sun(month))];                            %position of sun (km)

%----Radiation Network Matrix----%
Res_self_inv = ((1-emis_ir)./(A.*emis_ir)).^-1;                   
Res_other_inv = (repmat(transpose(A),1,ntot).*rad_vf);
Res_mat = Res_other_inv; 

for j=1:ntot 
    Res_mat(j,j)=-Res_self_inv(j)-sum(Res_other_inv(j,:));
end

%----Set Inital Time----%
temp_old=reshape(temp_init,1,1,ntot);
t_old=0;
delta_t_old=delta_t;
k_old=0;
N = t_final./delta_t;

%---Initialize arrays, to be filled in loop...----%
t = zeros(floor(N+1),1,1);                            
r_sc = zeros(floor(N+1),3,1);                         
v_sc = zeros(floor(N+1),3,1);                         
S = zeros(floor(N+1),3,1);                             
n = zeros(floor(N+1),3,ntot);                         
qext_arr = zeros(floor(N+1),3,ntot);                  
Qext = zeros(floor(N+1),1,ntot);
Qspace = zeros(floor(N+1),1,ntot);
Qcond = zeros(floor(N+1),1,ntot);
Qrad = zeros(floor(N+1),1,ntot);
Qint = zeros(floor(N+1),1,ntot);
temp = zeros(floor(N+1),1,ntot);
delta_t_lim = zeros(floor(N+1),1,ntot);
delta_t_lim_min = zeros(floor(N+1),1,1);
tcheck_sum = zeros(floor(N+1),1,1);
delta_t_new = zeros(floor(N+1),1,ntot);

nu1=0;                                                              
t1=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=======================Temperature Calculation Loop======================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while t_old < t_final
    
k=k_old+1;

%----Define Position and Velocity----%
t(k,:) = t_old+delta_t_old;                                       %current time of calculation(s)
[r_sc(k,:),v_sc(k,:)] = SCposition(acend_node,i,w,a,e,t(k,:));    %SC postion and velcocity

r_mag(k,:)=sqrt(r_sc(1)^2+r_sc(2)^2+r_sc(3)^2);

perp_2_orbit_plane = cross(r_sc(k,:),v_sc(k,:));                  %Vector perpendicular to orbit plane

S(k,:)=r_sun-r_sc(k,:);                                           %line of sight vector (SC and sun)


[EH1,EH2,nu2]=kepler(t1,nu1,t_old,a,e);                           %Resulting true anomaly from current timebnu
t1=t_old;
nu1=nu2;


[n1,n2,n3,n4,n5,n6]=normal(acend_node,i,w,nu2);                   %populate normal vector matrix
n(k,:,1)=[n1(1) n1(2) n1(3)];
n(k,:,2)=[n2(1) n2(2) n2(3)];
n(k,:,3)=[n3(1) n3(2) n3(3)];
n(k,:,4)=[n4(1) n4(2) n4(3)];
n(k,:,5)=[n5(1) n5(2) n5(3)];
n(k,:,6)=[n6(1) n6(2) n6(3)];


%--------------------------CALCULATION THROUGH NODES----------------------%
for j=1:ntot 

%----Calculate Heat Flow----%


%environmental radiation flux [qs,qa,qe] (W/m~2)
qext_arr(k,:,j) = ext_heat_flux_loop(r_sc(k,:),r_sun, n(k,:,j),Gs,AF,abso_sol(j),emis_ir(j),refl_cb(j),r_mag(k,:));
Qext_arr(k,:,j) = A(j).*qext_arr(k,:,j);                            

%environmental radiation input on the node (W)
Qext(k,:,j)= A(j).*sum(qext_arr(k,:,j),2);                          %Total external enviromental heat rate

%radiation to space
Qspace(k,:,j)=A(j).*emis_ir(j).*sigma*(0-temp_old(:,:,j).^4);       %Radiation to space 

%array of the node temp, differences (linear and 4th power)
temp_dif_mat = reshape((temp_old-temp_old(:,:,j)),1,ntot);
temp4_dif_mat = reshape((temp_old.^4-(temp_old(:,:,j)).^4),1,ntot);

%Conduction Heat flow
Qcond(k,:,j) = sum(conductance(j,:).*temp_dif_mat);                         

%%% Radiation Network %%%
E_mat = transpose(sigma.*(reshape(temp_old,1,ntot)).^4.*Res_self_inv);
J_mat = Res_mat\(-E_mat);
Q_rad(k,:,j) = sum(transpose(J_mat-J_mat(j,:)).*Res_other_inv(j,:));  %Radiative heat 
Qint(k,:,j)= inthl(j);                                                      %Internal Heat Load


%%%Check time-step for stability%%%
tcheck_cond = 1./C(j).*conductance(j,:);
ind = find(temp_dif_mat==0);
tcheck_rad = 1./C(j).*(transpose(J_mat-J_mat(j,:)).*Res_other_inv(j,:))./(temp_dif_mat);
tcheck_rad(ind) = 0;
tcheck_ext = 1./C(j).*Qext(k,:,j)./temp_old(:,:,j);
tcheck_space = 1./C(j).*Qspace(k,:,j)./(-temp_old(:,:,j));
tcheck_sum(k,:,j) = sum(tcheck_cond)+sum(abs(tcheck_rad))+sum(tcheck_ext)+sum(tcheck_space);
delta_t_lim(k,:,j) = 1./tcheck_sum(k,:,j);
end

delta_t_lim_min(k) = min(delta_t_lim(k,: ,:));          %time-step limit
    if delta_t > delta_t_lim_min(k)                     %if time-step too large
            delta_t_new = 0.9.*delta_t_lim_min(k);
    else
            delta_t_new = delta_t;
    end
    
%New temperatures
temp(k,:,:) = temp_old(:,:,:)+delta_t./(reshape(C,1,1,ntot)).*(Qext(k,:,:)+Qspace(k,:,:)+Qint(k,:,:)+Q_rad(k,:,:));

%Advance for next time loop
t_old=t(k,:);
delta_t_old = delta_t_new;
temp_old = temp(k,:,:);
k_old=k;

end

figure(1)
plot3(r_sc(:,1),r_sc(:,2),r_sc(:,3))
hold on
[x,y,z] = sphere();
surf( Re*x, Re*y, Re*z)
colormap(gray)
shading interp
hold on
sa = 2; st = 2.5;
zero_axis = [0 0 0];
x_axis = sa*Re*[1 0 0];
y_axis = sa*Re*[0 1 0];
z_axis = sa*Re*[0 0 1];
hold on
plot3(x_axis,zero_axis,zero_axis,'--k')
plot3(zero_axis,y_axis,zero_axis,'--k')
plot3(zero_axis,zero_axis,z_axis,'--k')
hold on
text(st*Re,0,0,'X'), text(0,st*Re,0,'Y'), text(0,0,st*Re,'Z')
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
title('Spacecraft Orbit')
axis square
grid on

figure(2)
subplot(2,3,1)
plot(t,Qext_arr(:,1,1))
hold on
plot(t,Qext_arr(:,2,1))
hold on
plot(t,Qext_arr(:,3,1))
hold on
xlabel('time (s)')
ylabel('W')
title('Face 1')
hold on
subplot(2,3,2)
plot(t,Qext_arr(:,1,2))
hold on
plot(t,Qext_arr(:,2,2))
hold on
plot(t,Qext_arr(:,3,2))
hold on
xlabel('time (s)')
ylabel('W')
title('Face 2')
hold on
subplot(2,3,3)
plot(t,Qext_arr(:,1,3))
hold on
plot(t,Qext_arr(:,2,3))
hold on
plot(t,Qext_arr(:,3,3))
hold on
xlabel('time (s)')
ylabel('W')
title('Face 3')
legend('Direct Solar','Albedo','Earth radiation')
hold on
subplot(2,3,4)
plot(t,Qext_arr(:,1,4))
hold on
plot(t,Qext_arr(:,2,4))
hold on
plot(t,Qext_arr(:,3,4))
hold on
xlabel('time (s)')
ylabel('W')
title('Face 4')
hold on
subplot(2,3,5)
plot(t,Qext_arr(:,1,5))
hold on
plot(t,Qext_arr(:,2,5))
hold on
plot(t,Qext_arr(:,3,5))
hold on
xlabel('time (s)')
ylabel('W')
title('Face 5')
subplot(2,3,6)
hold on
plot(t,Qext_arr(:,1,6))
hold on
plot(t,Qext_arr(:,2,6))
hold on
plot(t,Qext_arr(:,3,6))
hold on
xlabel('time (s)')
ylabel('W')
title('Face 6')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for i=1:6
    
%subplot(2,3,i)
%plot(t,Qext_arr(:,1,i))
%xlabel('time (s)')
%ylabel('W')
%title('Node')
%legend('Solar')

%end

figure(3)

plot(t,temp(:,:,1))
hold on
plot(t,temp(:,:,2))
hold on
plot(t,temp(:,:,3))
hold on
plot(t,temp(:,:,4))
hold on
plot(t,temp(:,:,5))
hold on
plot(t,temp(:,:,6))
xlabel('time (s)')
ylabel('T (K)')
title('Face Temperatures')
legend('Face 1','Face 2','Face 3','Face 4','Face 5','Face 6')