function [curve, tang, norm, bino]=curve_vdefortang_curvatureof1_orthogonalityerror_function(s,s0)

kap=1;

% T_solcoeff: matrix (array of coefficients of F_Tm which are the solutions to the tangent vector (T) VDE )(EQ. 4.5)
T_solcoeff=[[1,0,0, -(kap^2)/s0,(2*(kap^2))/s0^2, (kap^4/s0)*(1+(1/(s0^2))-(6/((kap*s0)^2)))];[0,1,0, -kap^2*(1+(1/s0^2)),(3*kap^2)/(s0^3),kap^4*(1+(3/(s0^2))+(2/(s0^4)))-((11*kap^2)/(s0^4))];[0,0,1,-2/s0,2*kap^2*((2/(kap*s0)^2)-1-(1/s0^2)),((2*kap^2)/(s0^3))*(6-(4/(kap^2))+(1/s0^2)-(2/((kap*s0)^2)))]];

% Function coefficients of tangent vector
F_T0= T_solcoeff(1,1) + T_solcoeff(1,2)*(s-s0)+ (1/factorial(2))*T_solcoeff(1,3)*(s-s0).^2 + (1/factorial(3))*T_solcoeff(1,4)*(s-s0).^3 + (1/factorial(4))*T_solcoeff(1,5)*(s-s0).^4 + (1/factorial(5))*T_solcoeff(1,6)*(s-s0).^5;
F_T1= T_solcoeff(2,1) + T_solcoeff(2,2)*(s-s0)+ (1/factorial(2))*T_solcoeff(2,3)*(s-s0).^2 + (1/factorial(3))*T_solcoeff(2,4)*(s-s0).^3 + (1/factorial(4))*T_solcoeff(2,5)*(s-s0).^4 + (1/factorial(5))*T_solcoeff(2,6)*(s-s0).^5;
F_T2= T_solcoeff(3,1) + T_solcoeff(3,2)*(s-s0)+ (1/factorial(2))*T_solcoeff(3,3)*(s-s0).^2 + (1/factorial(3))*T_solcoeff(3,4)*(s-s0).^3 + (1/factorial(4))*T_solcoeff(3,5)*(s-s0).^4 + (1/factorial(5))*T_solcoeff(3,6)*(s-s0).^5;

% Function coefficients of normal vector
F_N0= (T_solcoeff(1,2)+ T_solcoeff(1,3)*(s-s0) + (1/factorial(2))*T_solcoeff(1,4)*(s-s0).^2 + (1/factorial(3))*T_solcoeff(1,5)*(s-s0).^3 + (1/factorial(4))*T_solcoeff(1,6)*(s-s0).^4)/kap;
F_N1= (T_solcoeff(2,2)+ T_solcoeff(2,3)*(s-s0) + (1/factorial(2))*T_solcoeff(2,4)*(s-s0).^2 + (1/factorial(3))*T_solcoeff(2,5)*(s-s0).^3 + (1/factorial(4))*T_solcoeff(2,6)*(s-s0).^4)/kap;
F_N2= (T_solcoeff(3,2)+ T_solcoeff(3,3)*(s-s0) + (1/factorial(2))*T_solcoeff(3,4)*(s-s0).^2 + (1/factorial(3))*T_solcoeff(3,5)*(s-s0).^3 + (1/factorial(4))*T_solcoeff(3,6)*(s-s0).^4)/kap;

% Function coefficients of curve vector
F_curve0= T_solcoeff(1,1)*(s-s0) + (1/factorial(2))*T_solcoeff(1,2)*(s-s0).^2+ (1/factorial(3))*T_solcoeff(1,3)*(s-s0).^3 + (1/factorial(4))*T_solcoeff(1,4)*(s-s0).^4 + (1/factorial(5))*T_solcoeff(1,5)*(s-s0).^5 + (1/factorial(6))*T_solcoeff(1,6)*(s-s0).^6;
F_curve1= T_solcoeff(2,1)*(s-s0) + (1/factorial(2))*T_solcoeff(2,2)*(s-s0).^2+ (1/factorial(3))*T_solcoeff(2,3)*(s-s0).^3 + (1/factorial(4))*T_solcoeff(2,4)*(s-s0).^4 + (1/factorial(5))*T_solcoeff(2,5)*(s-s0).^5 + (1/factorial(6))*T_solcoeff(2,6)*(s-s0).^6;
F_curve2= T_solcoeff(3,1)*(s-s0) + (1/factorial(2))*T_solcoeff(3,2)*(s-s0).^2+ (1/factorial(3))*T_solcoeff(3,3)*(s-s0).^3 + (1/factorial(4))*T_solcoeff(3,4)*(s-s0).^4 + (1/factorial(5))*T_solcoeff(3,5)*(s-s0).^5 + (1/factorial(6))*T_solcoeff(3,6)*(s-s0).^6;

% Tangent vector
tang(1,:)=F_T0;
tang(2,:)=F_T1;
tang(3,:)=F_T2;

% Principal normal vector
norm(1,:)=F_N0;
norm(2,:)=F_N1;
norm(3,:)=F_N2;

% Binormal vector
bino(1,:)=tang(2,:).*norm(3,:)-tang(3,:).*norm(2,:);
bino(2,:)=tang(3,:).*norm(1,:)-tang(1,:).*norm(3,:);
bino(3,:)=tang(1,:).*norm(2,:)-tang(2,:).*norm(1,:);

% Curve
curve(1,:)=F_curve0;
curve(2,:)=F_curve1;
curve(3,:)=F_curve2;
end