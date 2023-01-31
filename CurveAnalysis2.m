% In this script, we look at the curves obtained from solving the VDE for 
% the tangent vector but the constant vectors a0, a1, a2 are assumed to be
% orthogonal.
close all hidden % removes all current figures 
clear % clears any variables in your workspace

% GENERAL PARAMETERS/INFORMATION
% arc length domain
s0=5; % initial arclength (VARY)
s_int=0.9; % size of interval for s (VARY)
sf=s0+s_int; % final arclength value (DO NOT CHANGE)
ds=0.00002; % step size (VARY)
s=s0:ds:sf; % arc length array (DO NOT CHANGE)
kap=2; % curvature (VARY)

%% Looking at orthogonality issue for curvature=1, s0=any. Curves are based on work from third draft

% Obtaining data for curves (s0=1 & curvature=1)
[curve, tang, norm, bino]=curve_vdefortang_curvatureof1_function(s,s0);
[curve_oe, tang_oe, norm_oe, bino_oe]=curve_vdefortang_curvatureof1_orthogonalityerror_function(s,s0);

% Calculating magnitude of tangent and normal vectors
magn(1,:)=sqrt(tang(1,:).^2 + tang(2,:).^2 + tang(3,:).^2);
magn(2,:)=sqrt(norm(1,:).^2 + norm(2,:).^2 + norm(3,:).^2);
magn(3,:)=sqrt(bino(1,:).^2 + bino(2,:).^2 + bino(3,:).^2);
magn(4,:)=sqrt(tang_oe(1,:).^2 + tang_oe(2,:).^2 + tang_oe(3,:).^2);
magn(5,:)=sqrt(norm_oe(1,:).^2 + norm_oe(2,:).^2 + norm_oe(3,:).^2);
magn(6,:)=sqrt(bino_oe(1,:).^2 + bino_oe(2,:).^2 + bino_oe(3,:).^2);

% Calculating dot product of tangent and normal vectors
dot_tn(1,:)= tang(1,:).*norm(1,:) + tang(2,:).*norm(2,:) + tang(3,:).*norm(3,:);
dot_tn(2,:)= tang_oe(1,:).*norm_oe(1,:) + tang_oe(2,:).*norm_oe(2,:) + tang_oe(3,:).*norm_oe(3,:);

%--------------------------------------------------------------------------------
%PLOTTING
%--------------------------------------------------------------------------------

% Plotting curves
figure(1)
camproj('orthographic')
hold off
plot3(curve(1,:), curve(2,:), curve(3,:), "r",'linewidth',2)
axis equal
grid on
hold on
plot3(curve(1,1), curve(2,1), curve(3,1), '-o','Color','b','MarkerSize',8,'MarkerFaceColor','#D9FFFF')
plot3(curve_oe(1,:), curve_oe(2,:), curve_oe(3,:), "b",'linewidth',2)
xlabel("x")
ylabel("y")
zlabel("z")

% Plotting magnitude of tangent vector
figure(2)
plot(s, magn(1,:),"r")
hold on
plot(s, magn(4,:),"b")
hold off

% Plotting magnitude of principal normal vector
figure(3)
plot(s, magn(2,:),"r")
hold on
plot(s, magn(5,:),"b")
hold off

% Plotting magnitude of binormal vector
figure(4)
plot(s, magn(3,:),"r")
hold on
plot(s, magn(6,:),"b")
hold off


% Plotting orthogonality of tangent and principal normal vectors
figure(5)
plot(s, dot_tn(1,:),"r")
hold on
plot(s, dot_tn(2,:),"b")
hold off

%% Looking at orthogonality issue for curvature=any, s0=1. Curves are based on work from second draft

% Obtaining data for curves (s0=1 & curvature=1)
[curve1, tang1, norm1, bino1]=curve_vdefornorm_nonexpanded_function(s,s0);
[curve1_oe, tang1_oe, norm1_oe, bino1_oe]=curve_vdefornorm_nonexpanded_orthogonalityerror_function(s,s0);


% Calculating magnitude of tangent and normal vectors
magn1(1,:)=sqrt(tang1(1,:).^2 + tang1(2,:).^2 + tang1(3,:).^2);
magn1(2,:)=sqrt(norm1(1,:).^2 + norm1(2,:).^2 + norm1(3,:).^2);
magn1(3,:)=sqrt(bino1(1,:).^2 + bino1(2,:).^2 + bino1(3,:).^2);
magn1(4,:)=sqrt(tang1_oe(1,:).^2 + tang1_oe(2,:).^2 + tang1_oe(3,:).^2);
magn1(5,:)=sqrt(norm1_oe(1,:).^2 + norm1_oe(2,:).^2 + norm1_oe(3,:).^2);
magn1(6,:)=sqrt(bino1_oe(1,:).^2 + bino1_oe(2,:).^2 + bino1_oe(3,:).^2);

% Calculating dot product of tangent and normal vectors
dot_tn1(1,:)= tang1(1,:).*norm1(1,:) + tang1(2,:).*norm1(2,:) + tang1(3,:).*norm1(3,:);
dot_tn1(2,:)= tang1_oe(1,:).*norm1_oe(1,:) + tang1_oe(2,:).*norm1_oe(2,:) + tang1_oe(3,:).*norm1_oe(3,:);


%{

% Plotting curves
figure(6)
camproj('orthographic')
hold off
plot3(curve1(1,:), curve1(2,:), curve1(3,:), "r",'linewidth',2)
axis equal
grid on
hold on
plot3(curve1(1,1), curve1(2,1), curve1(3,1), '-o','Color','b','MarkerSize',8,'MarkerFaceColor','#D9FFFF')
plot3(curve1_oe(1,:), curve1_oe(2,:), curve1_oe(3,:), "b",'linewidth',2)
xlabel("x")
ylabel("y")
zlabel("z")

% Plotting magnitude of tangent vector
figure(7)
plot(s, magn1(1,:),"r")
hold on
plot(s, magn1(4,:),"b")
hold off

% Plotting magnitude of principal normal vector
figure(8)
plot(s, magn1(2,:),"r")
hold on
plot(s, magn1(5,:),"b")
hold off

% Plotting magnitude of binormal vector
figure(9)
plot(s, magn1(3,:),"r")
hold on
plot(s, magn1(6,:),"b")
hold off

% Plotting orthogonality of tangent and principal normal vectors
figure(10)
plot(s, dot_tn1(1,:),"r")
hold on
plot(s, dot_tn1(2,:),"b")
hold off
%}