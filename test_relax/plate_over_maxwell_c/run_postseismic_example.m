%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script to illustrate use of PLate_over_Maxwell codes

%fault geometry  -- type 'help Plate_over_Maxwell_Layer_over_Halfspace_Displacements' for
%more info
%m=[length,width,depth of bottom edge,dip,strike,east*,north*,strike-slip,dip-slip,opening]
% *position of midpoint of bottom edge
m=[200 20 20 90 0 0 0 1 0 0];


%coordinates of points on surface
x=linspace(-300,300,24);
y=linspace(-300,300,24);
[X,Y]=meshgrid(x,y);
xloc=[X(:)';Y(:)'];

%material properties
H1 = 20;  %elastic thickness (km)
H2 = 40;  %depth to bottom of viscoelastic layer (km)
tR1 = 25;  %relaxation time of visco layer (years)
tR2 = 25;   %relaxation time of halfspace (years)
nu=0.25;   %Poisson's ratio


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute displacements and velocities at different times

time=1;  %time since EQ (years)
%U is a 3xn matrix of displacements (first row east, second row north, third row up) 
%note: cumulative postseismic displacement (no coseismic)
U_1 = Plate_over_Maxwell_Layer_over_Halfspace_Displacements(m,xloc,H1,H2,nu,time,tR1,tR2);
V_1 = Plate_over_Maxwell_Layer_over_Halfspace_Velocities(m,xloc,H1,H2,nu,time,tR1,tR2);


time=5;  %time since EQ (years)
U_5 = Plate_over_Maxwell_Layer_over_Halfspace_Displacements(m,xloc,H1,H2,nu,time,tR1,tR2);
V_5 = Plate_over_Maxwell_Layer_over_Halfspace_Velocities(m,xloc,H1,H2,nu,time,tR1,tR2);


time=10;  %time since EQ (years)
U_10 = Plate_over_Maxwell_Layer_over_Halfspace_Displacements(m,xloc,H1,H2,nu,time,tR1,tR2);
V_10 = Plate_over_Maxwell_Layer_over_Halfspace_Velocities(m,xloc,H1,H2,nu,time,tR1,tR2);


time=50;   %time since EQ (years)
U_50 = Plate_over_Maxwell_Layer_over_Halfspace_Displacements(m,xloc,H1,H2,nu,time,tR1,tR2);
V_50 = Plate_over_Maxwell_Layer_over_Halfspace_Velocities(m,xloc,H1,H2,nu,time,tR1,tR2);


time=100;   %time since EQ (years)
U_100 = Plate_over_Maxwell_Layer_over_Halfspace_Displacements(m,xloc,H1,H2,nu,time,tR1,tR2);
V_100 = Plate_over_Maxwell_Layer_over_Halfspace_Velocities(m,xloc,H1,H2,nu,time,tR1,tR2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot displacement or velocity vectors
figure
scale = 500;
quiver(xloc(1,:),xloc(2,:),scale*U_1(1,:),scale*U_1(2,:),0)
hold on
quiver(xloc(1,:),xloc(2,:),scale*U_5(1,:),scale*U_5(2,:),0)
quiver(xloc(1,:),xloc(2,:),scale*U_10(1,:),scale*U_10(2,:),0)
quiver(xloc(1,:),xloc(2,:),scale*U_50(1,:),scale*U_50(2,:),0)
quiver(xloc(1,:),xloc(2,:),scale*U_100(1,:),scale*U_100(2,:),0)
axis equal
title('cumulative postseismic displacements')
legend('t=1','t=5','t=10','t=50','t=100')

figure
scale = 10000;
quiver(xloc(1,:),xloc(2,:),scale*V_1(1,:),scale*V_1(2,:),0)
hold on
quiver(xloc(1,:),xloc(2,:),scale*V_5(1,:),scale*V_5(2,:),0)
quiver(xloc(1,:),xloc(2,:),scale*V_10(1,:),scale*V_10(2,:),0)
quiver(xloc(1,:),xloc(2,:),scale*V_50(1,:),scale*V_50(2,:),0)
quiver(xloc(1,:),xloc(2,:),scale*V_100(1,:),scale*V_100(2,:),0)
axis equal
title('postseismic velocities')
legend('t=1','t=5','t=10','t=50','t=100')
