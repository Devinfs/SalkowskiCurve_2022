close all hidden % removes all current figures 
clear % clears any variables in your workspace

% GENERAL PARAMETERS/INFORMATION
% arc length domain
s0=1; % initial arclength (DO NOT CHANGE)
s_int=3.9; % size of interval for s (VARY)
sf=s0+s_int; % final arclength value (DO NOT CHANGE)
ds=0.00002; % step size (VARY)
s=s0:ds:sf; % arc length array (DO NOT CHANGE)
kap=4; % curvature (VARY)

%% Looking at integration issue for curvature=1, s0=any. Curves are based on work from clarification draft

% Defining integration constants for work
integrationconstant1=[0,0,0];
integrationconstant2=[0,0,0];

% Obtaining data for curves 
[curve_ex, tang_ex, norm_ex, bino_ex]=curve_vdefornorm_expanded_function(s,kap);
[curve_nex, tang_nex, norm_nex, bino_nex]=curve_vdefornorm_nonexpanded_function(s,kap);
[curve_ex_ii, tang_ex_ii, norm_ex_ii, bino_ex_ii]=curve_vdefornorm_expanded_integrationissue_function(s,kap, integrationconstant1, integrationconstant2);
[curve_nex_ii, tang_nex_ii, norm_nex_ii, bino_nex_ii]=curve_vdefornorm_nonexpanded_integrationissue_function(s,kap, integrationconstant1, integrationconstant2 );


% Calculating magnitude of Frenet frame vectors
magn(1,:)=sqrt(tang_ex(1,:).^2 + tang_ex(2,:).^2 + tang_ex(3,:).^2);
magn(2,:)=sqrt(tang_nex(1,:).^2 + tang_nex(2,:).^2 + tang_nex(3,:).^2);
magn(3,:)=sqrt(tang_ex_ii(1,:).^2 + tang_ex_ii(2,:).^2 + tang_ex_ii(3,:).^2);
magn(4,:)=sqrt(tang_nex_ii(1,:).^2 + tang_nex_ii(2,:).^2 + tang_nex_ii(3,:).^2);
magn(5,:)=sqrt(norm_ex(1,:).^2 + norm_ex(2,:).^2 + norm_ex(3,:).^2);
magn(6,:)=sqrt(norm_nex(1,:).^2 + norm_nex(2,:).^2 + norm_nex(3,:).^2);
magn(7,:)=sqrt(norm_ex_ii(1,:).^2 + norm_ex_ii(2,:).^2 + norm_ex_ii(3,:).^2);
magn(8,:)=sqrt(norm_nex_ii(1,:).^2 + norm_nex_ii(2,:).^2 + norm_nex_ii(3,:).^2);
magn(9,:)=sqrt(bino_ex(1,:).^2 + bino_ex(2,:).^2 + bino_ex(3,:).^2);
magn(10,:)=sqrt(bino_nex(1,:).^2 + bino_nex(2,:).^2 + bino_nex(3,:).^2);
magn(11,:)=sqrt(bino_ex_ii(1,:).^2 + bino_ex_ii(2,:).^2 + bino_ex_ii(3,:).^2);
magn(12,:)=sqrt(bino_nex_ii(1,:).^2 + bino_nex_ii(2,:).^2 + bino_nex_ii(3,:).^2);

%--------------------------------------------------------------------------------
%PLOTTING
%--------------------------------------------------------------------------------


% Plotting magnitude of tangent vector
figure(2)
plot(s, magn(1,:),"r")
hold on
plot(s, magn(2,:),"b")
plot(s, magn(3,:),"g")
plot(s, magn(4,:),"k")
hold off


% Plotting magnitude of principal normal vector
figure(3)
plot(s, magn(5,:),"r")
hold on
plot(s, magn(6,:),"b")
plot(s, magn(7,:),"g")
plot(s, magn(8,:),"k")
hold off


% Plotting magnitude of binormal vector
figure(4)
plot(s, magn(9,:),"r")
hold on
plot(s, magn(10,:),"b")
plot(s, magn(11,:),"g")
plot(s, magn(12,:),"k")
hold off