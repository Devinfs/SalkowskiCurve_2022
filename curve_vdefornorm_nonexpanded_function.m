function [curve, tang, norm, bino]=curve_vdefornorm_nonexpanded_function(s, kap)


% initial s
s0=1;
% t interval
t=kap*log(s);

% B matrix (array of coefficients of F_Nm)(EQ. 4.5)
B=[[1,0,0, -1/(factorial(3)*kap),-1/(factorial(4)*(kap^2)),(1/factorial(5))*((-13/(kap^3))+(1/kap))];[0,1,0, -2/factorial(3),-5/(factorial(4)*kap),(1/factorial(5))*((-17/(kap^2))+4)];[0,0,1,2/(factorial(3)*kap),(1/factorial(4))*((2/(kap^2))-4),(1/factorial(5))*((2/(kap^3))-(18/kap))]];

% functions needed for normal vector
F_N0= B(1,1) + B(1,2)*t + + B(1,3)*t.^2 + B(1,4)*t.^3 + B(1,5)*t.^4 + B(1,6)*t.^5;
F_N1= B(2,1) + B(2,2)*t + + B(2,3)*t.^2 + B(2,4)*t.^3 + B(2,5)*t.^4 + B(2,6)*t.^5;
F_N2= B(3,1) + B(3,2)*t + + B(3,3)*t.^2 + B(3,4)*t.^3 + B(3,5)*t.^4 + B(3,6)*t.^5;

% PRINCIPAL NORMAL VECTOR (SAME FOR BOTH METHODS)
norm(1,:)=-F_N1 - ((1/(2*kap))*F_N2);
norm(2,:)=F_N0 - F_N2;
norm(3,:)=F_N1;


% The next two "%" relate to calculating C_mk coefficients from Section 6.1
% EN matrix (array with values needed to calculate Cmk values)
EN_nonexp=[[1,0,0,1,1,1];[0,0,0,1,1,1];[0,0,0,1,1,1];[0,0,0,1,1,1];[0,0,0,0,1,1];[0,0,0,0,0,1]];
EN_nonexp(:,:,2)=[[0,1,0,1,1,1];[0,1,0,1,1,1];[0,0,0,1,1,1];[0,0,0,1,1,1];[0,0,0,0,1,1];[0,0,0,0,0,1]];
EN_nonexp(:,:,3)=[[0,0,1,1,1,1];[0,0,1,1,1,1];[0,0,1,1,1,1];[0,0,0,1,1,1];[0,0,0,0,1,1];[0,0,0,0,0,1]];

% Cmk matrix (array of coefficients of F_Tm)
Cmk=zeros(3,6);
for j=0:5
    % Calculating coefficients of first function
    Cmk(1,1)=Cmk(1,1) + ((-1)^(j))*factorial(j)*(kap^(j))*EN_nonexp(1,j+1,1)*B(1,j+1);
    Cmk(1,2)=Cmk(1,2) + ((-1)^(j-1))*factorial(j)*(kap^(j-1))*EN_nonexp(2,j+1,1)*B(1,j+1);
    Cmk(1,3)=Cmk(1,3) + ((-1)^(j))*(factorial(j)/factorial(2))*(kap^(j-2))*EN_nonexp(3,j+1,1)*B(1,j+1);
    Cmk(1,4)=Cmk(1,4) + ((-1)^(j-3))*(factorial(j)/factorial(3))*(kap^(j-3))*EN_nonexp(4,j+1,1)*B(1,j+1);
    Cmk(1,5)=Cmk(1,5) + ((-1)^(j-4))*(factorial(j)/factorial(4))*(kap^(j-4))*EN_nonexp(5,j+1,1)*B(1,j+1);
    Cmk(1,6)=Cmk(1,6) + ((-1)^(j-5))*(factorial(j)/factorial(5))*(kap^(j-5))*EN_nonexp(6,j+1,1)*B(1,j+1);
    % Calculating coefficients of second function
    Cmk(2,1)=Cmk(2,1) + ((-1)^(j))*factorial(j)*(kap^(j))*EN_nonexp(1,j+1,2)*B(2,j+1);
    Cmk(2,2)=Cmk(2,2) + ((-1)^(j-1))*factorial(j)*(kap^(j-1))*EN_nonexp(2,j+1,2)*B(2,j+1);
    Cmk(2,3)=Cmk(2,3) + ((-1)^(j))*(factorial(j)/factorial(2))*(kap^(j-2))*EN_nonexp(3,j+1,2)*B(2,j+1);
    Cmk(2,4)=Cmk(2,4) + ((-1)^(j-3))*(factorial(j)/factorial(3))*(kap^(j-3))*EN_nonexp(4,j+1,2)*B(2,j+1);
    Cmk(2,5)=Cmk(2,5) + ((-1)^(j-4))*(factorial(j)/factorial(4))*(kap^(j-4))*EN_nonexp(5,j+1,2)*B(2,j+1);
    Cmk(2,6)=Cmk(2,6) + ((-1)^(j-5))*(factorial(j)/factorial(5))*(kap^(j-5))*EN_nonexp(6,j+1,2)*B(2,j+1);
    % Calculating coefficients of third function
    Cmk(3,1)=Cmk(3,1) + ((-1)^(j))*factorial(j)*(kap^(j))*EN_nonexp(1,j+1,3)*B(3,j+1);
    Cmk(3,2)=Cmk(3,2) + ((-1)^(j-1))*factorial(j)*(kap^(j-1))*EN_nonexp(2,j+1,3)*B(3,j+1);
    Cmk(3,3)=Cmk(3,3) + ((-1)^(j))*(factorial(j)/factorial(2))*(kap^(j-2))*EN_nonexp(3,j+1,3)*B(3,j+1);
    Cmk(3,4)=Cmk(3,4) + ((-1)^(j-3))*(factorial(j)/factorial(3))*(kap^(j-3))*EN_nonexp(4,j+1,3)*B(3,j+1);
    Cmk(3,5)=Cmk(3,5) + ((-1)^(j-4))*(factorial(j)/factorial(4))*(kap^(j-4))*EN_nonexp(5,j+1,3)*B(3,j+1);
    Cmk(3,6)=Cmk(3,6) + ((-1)^(j-5))*(factorial(j)/factorial(5))*(kap^(j-5))*EN_nonexp(6,j+1,3)*B(3,j+1);
end

% functions needed for tangent vector (non-expanded method)
F_T0= kap*exp(t./kap).*(Cmk(1,1) + Cmk(1,2)*t + + Cmk(1,3)*t.^2 + Cmk(1,4)*t.^3 + Cmk(1,5)*t.^4 + Cmk(1,6)*t.^5) - kap*Cmk(1,1);
F_T1= kap*exp(t./kap).*(Cmk(2,1) + Cmk(2,2)*t + + Cmk(2,3)*t.^2 + Cmk(2,4)*t.^3 + Cmk(2,5)*t.^4 + Cmk(2,6)*t.^5) - kap*Cmk(2,1);
F_T2= kap*exp(t./kap).*(Cmk(3,1) + Cmk(3,2)*t + + Cmk(3,3)*t.^2 + Cmk(3,4)*t.^3 + Cmk(3,5)*t.^4 + Cmk(3,6)*t.^5) - kap*Cmk(3,1);

% TANGENT VECTOR (NON-EXPANDED METHOD)
tang(1,:)=1-F_T1 - (1/(2*kap))*F_T2;
tang(2,:)= F_T0 - F_T2;
tang(3,:)=F_T1;


% BINORMAL VECTOR (NON-EXPANDED METHOD)
bino(1,:)=tang(2,:).*norm(3,:)-tang(3,:).*norm(2,:);
bino(2,:)=tang(3,:).*norm(1,:)-tang(1,:).*norm(3,:);
bino(3,:)=tang(1,:).*norm(2,:)-tang(2,:).*norm(1,:);


% The next two "%" relate to calculating D_ml coefficients from Appendix B
% ET matrix (array with values needed to calculate Dqk values)
ET_nonexp=[[1/2, 1/4, 1/4, 3/8, 3/4, 15/8];[0, 1/2, 1/2, 3/4, 3/2, 15/4];[0, 0, 1/2, 3/4, 3/2, 15/4];[0, 0, 0, 1/2, 1, 5/2];[0, 0, 0, 0, 1/2, 5/4];[0, 0, 0, 0, 0, 1/2]];

% D_ml matrix (array of coefficients of curve function terms)
D_ml=zeros(3,6);
for j=0:5
    % Calculating coefficients of first function
    D_ml(1,1)=D_ml(1,1)+ ((-1)^(j))*(kap^(j))*ET_nonexp(1,j+1)*Cmk(1,j+1);
    D_ml(1,2)=D_ml(1,2)+ ((-1)^(j-1))*(kap^(j))*ET_nonexp(2,j+1)*Cmk(1,j+1);
    D_ml(1,3)=D_ml(1,3)+ ((-1)^(j-2))*(kap^(j))*ET_nonexp(3,j+1)*Cmk(1,j+1);
    D_ml(1,4)=D_ml(1,4)+ ((-1)^(j-3))*(kap^(j))*ET_nonexp(4,j+1)*Cmk(1,j+1);
    D_ml(1,5)=D_ml(1,5)+ ((-1)^(j-4))*(kap^(j))*ET_nonexp(5,j+1)*Cmk(1,j+1);
    D_ml(1,6)=D_ml(1,6)+ ((-1)^(j-5))*(kap^(j))*ET_nonexp(6,j+1)*Cmk(1,j+1);
    % Calculating coefficients of second function
    D_ml(2,1)=D_ml(2,1)+ ((-1)^(j))*(kap^(j))*ET_nonexp(1,j+1)*Cmk(2,j+1);
    D_ml(2,2)=D_ml(2,2)+ ((-1)^(j-1))*(kap^(j))*ET_nonexp(2,j+1)*Cmk(2,j+1);
    D_ml(2,3)=D_ml(2,3)+ ((-1)^(j-2))*(kap^(j))*ET_nonexp(3,j+1)*Cmk(2,j+1);
    D_ml(2,4)=D_ml(2,4)+ ((-1)^(j-3))*(kap^(j))*ET_nonexp(4,j+1)*Cmk(2,j+1);
    D_ml(2,5)=D_ml(2,5)+ ((-1)^(j-4))*(kap^(j))*ET_nonexp(5,j+1)*Cmk(2,j+1);
    D_ml(2,6)=D_ml(2,6)+ ((-1)^(j-5))*(kap^(j))*ET_nonexp(6,j+1)*Cmk(2,j+1);
    % Calculating coefficients of third function
    D_ml(3,1)=D_ml(3,1)+ ((-1)^(j))*(kap^(j))*ET_nonexp(1,j+1)*Cmk(3,j+1);
    D_ml(3,2)=D_ml(3,2)+ ((-1)^(j-1))*(kap^(j))*ET_nonexp(2,j+1)*Cmk(3,j+1);
    D_ml(3,3)=D_ml(3,3)+ ((-1)^(j-2))*(kap^(j))*ET_nonexp(3,j+1)*Cmk(3,j+1);
    D_ml(3,4)=D_ml(3,4)+ ((-1)^(j-3))*(kap^(j))*ET_nonexp(4,j+1)*Cmk(3,j+1);
    D_ml(3,5)=D_ml(3,5)+ ((-1)^(j-4))*(kap^(j))*ET_nonexp(5,j+1)*Cmk(3,j+1);
    D_ml(3,6)=D_ml(3,6)+ ((-1)^(j-5))*(kap^(j))*ET_nonexp(6,j+1)*Cmk(3,j+1);
end


% The next two "%" relate to the final functions needed for curve (EQ. 9.1)
% integration constants of main functions
F_curve0_cons=kap*(s0.^2).*(D_ml(1,1)+ D_ml(1,2)*log(s0)+ D_ml(1,3)*(log(s0).^2) + D_ml(1,4)*(log(s0).^3) + D_ml(1,5)*(log(s0).^4) + + D_ml(1,6)*(log(s0).^5));
F_curve1_cons=kap*(s0.^2).*(D_ml(2,1)+ D_ml(2,2)*log(s0)+ D_ml(2,3)*(log(s0).^2) + D_ml(2,4)*(log(s0).^3) + D_ml(2,5)*(log(s0).^4) + + D_ml(2,6)*(log(s0).^5));
F_curve2_cons=kap*(s0.^2).*(D_ml(3,1)+ D_ml(3,2)*log(s0)+ D_ml(3,3)*(log(s0).^2) + D_ml(3,4)*(log(s0).^3) + D_ml(3,5)*(log(s0).^4) + + D_ml(3,6)*(log(s0).^5));

% main functions
F_curve0= kap*(s.^2).*(D_ml(1,1)+ D_ml(1,2)*log(s)+ D_ml(1,3)*(log(s).^2) + D_ml(1,4)*(log(s).^3) + D_ml(1,5)*(log(s).^4) + + D_ml(1,6)*(log(s).^5)) - F_curve0_cons - kap*Cmk(1,1)*(s-s0);
F_curve1= kap*(s.^2).*(D_ml(2,1)+ D_ml(2,2)*log(s)+ D_ml(2,3)*(log(s).^2) + D_ml(2,4)*(log(s).^3) + D_ml(2,5)*(log(s).^4) + + D_ml(2,6)*(log(s).^5)) - F_curve1_cons - kap*Cmk(2,1)*(s-s0);
F_curve2= kap*(s.^2).*(D_ml(3,1)+ D_ml(3,2)*log(s)+ D_ml(3,3)*((log(s).^2)) + D_ml(3,4)*(log(s).^3) + D_ml(3,5)*(log(s).^4) + + D_ml(3,6)*(log(s).^5)) - F_curve2_cons - kap*Cmk(3,1)*(s-s0);

% NON-EXPANDED METHOD
curve(1,:)= s-s0-(1/(2*kap))*F_curve2 - F_curve1; %COEFFICIENT OF VECTOR A
curve(2,:)=F_curve0-F_curve2;    %COEFFICIENT OF VECTOR B
curve(3,:)=-F_curve1; %COEFFICIENT OF VECTOR C

end