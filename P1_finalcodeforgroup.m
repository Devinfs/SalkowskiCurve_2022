% Code needed to plot calculate coefficients of functions to plot 
% curves from expanded and non-expanded methods.
clear

% GENERAL PARAMETERS/INFORMATION
% arc length domain
s0=1; % initial arclength (do not change)
sf=2.0; % final arclength value (can be changed)
s=s0:0.0001:sf;

% value of curvature
kap=2;

% t interval
t=kap*log(s);


% tangent vector coefficient functions (F_Tm)
% normal vector coefficient functions (F_Nm)

% B matrix (array of coefficients of F_Nm)
B=[[1,0,0, -1/(factorial(3)*kap),-1/(factorial(4)*(kap^2)),(1/factorial(5))*((-13/(kap^3))+(1/kap))];[0,1,0, -2/factorial(3),-5/(factorial(4)*kap),(1/factorial(5))*((-17/(kap^2))+4)];[0,0,1,2/(factorial(3)*kap),(1/factorial(4))*((2/(kap^2))-4),(1/factorial(5))*((2/(kap^3))-(18/kap))]];


% EXPANDED METHOD (expanding the exponential term in integral to obtain T)

% Cp matrix (array of coeffecients of F_Tm)
Cp= [[1, 1/(2*kap), 1/(factorial(3)*kap^2), (1/(factorial(4)))*((1/kap^3)-(1/kap)), (1/factorial(5))*((1/kap^4)-(9/kap^2)),0,0,0,0,0]; [0, 1/2, 1/(3*kap), ((1/(8*kap^2))-(1/12)),((1/(30*kap^3))-(13/(120*kap^2))),0,0,0,0,0];[0,0,1/3,1/(3*kap),((-1/30)+(11/(60*kap^2))), ((1/(45*kap^3))-(19/(360*kap)))],0,0,0,0];

% ET matrix (array with values needed to calculate Dp values)
ET_exp=[[1, 1, 1, 1, 1, 1]; [1, 1, 1, 1, 1, 1]; [0, 1, 1, 1, 1, 1]; [0, 0, 1, 1, 1, 1]; [0, 0, 0, 1, 1, 1]; [0, 0, 0, 0, 1, 1]; [0, 0, 0, 0, 0, 1]];

% Dp matrix (array of coefficients of main functions)
Dp=zeros(3,7);
for j=1:6
    % Calculating coefficients of first function
    Dp(1,1)=Dp(1,1)+ ((-1)^(j))*(kap^(j))*factorial(j)*ET_exp(1,j)*Cp(1,j);
    Dp(1,2)=Dp(1,2)+ ((-1)^(j-1))*(kap^(j))*factorial(j)*ET_exp(2,j)*Cp(1,j);
    Dp(1,3)=Dp(1,3)+ ((-1)^(j-2))*(kap^(j))*(factorial(j)/factorial(2))*ET_exp(3,j)*Cp(1,j);
    Dp(1,4)=Dp(1,4)+ ((-1)^(j-3))*(kap^(j))*(factorial(j)/factorial(3))*ET_exp(4,j)*Cp(1,j);
    Dp(1,5)=Dp(1,5)+ ((-1)^(j-4))*(kap^(j))*(factorial(j)/factorial(4))*ET_exp(5,j)*Cp(1,j);
    Dp(1,6)=Dp(1,6)+ ((-1)^(j-5))*(kap^(j))*(factorial(j)/factorial(5))*ET_exp(6,j)*Cp(1,j);
    Dp(1,7)=Dp(1,7)+ ((-1)^(j-6))*(kap^(j))*(factorial(j)/factorial(6))*ET_exp(7,j)*Cp(1,j);
    % Calculating coefficients of second function
    Dp(2,1)=Dp(2,1)+ ((-1)^(j))*(kap^(j))*factorial(j)*ET_exp(1,j)*Cp(2,j);
    Dp(2,2)=Dp(2,2)+ ((-1)^(j-1))*(kap^(j))*factorial(j)*ET_exp(2,j)*Cp(2,j);
    Dp(2,3)=Dp(2,3)+ ((-1)^(j-2))*(kap^(j))*(factorial(j)/factorial(2))*ET_exp(3,j)*Cp(2,j);
    Dp(2,4)=Dp(2,4)+ ((-1)^(j-3))*(kap^(j))*(factorial(j)/factorial(3))*ET_exp(4,j)*Cp(2,j);
    Dp(2,5)=Dp(2,5)+ ((-1)^(j-4))*(kap^(j))*(factorial(j)/factorial(4))*ET_exp(5,j)*Cp(2,j);
    Dp(2,6)=Dp(2,6)+ ((-1)^(j-5))*(kap^(j))*(factorial(j)/factorial(5))*ET_exp(6,j)*Cp(2,j);
    Dp(2,7)=Dp(2,7)+ ((-1)^(j-6))*(kap^(j))*(factorial(j)/factorial(6))*ET_exp(7,j)*Cp(2,j);
    % Calculating coefficients of third function
    Dp(3,1)=Dp(3,1)+ ((-1)^(j))*(kap^(j))*factorial(j)*ET_exp(1,j)*Cp(3,j);
    Dp(3,2)=Dp(3,2)+ ((-1)^(j-1))*(kap^(j))*factorial(j)*ET_exp(2,j)*Cp(3,j);
    Dp(3,3)=Dp(3,3)+ ((-1)^(j-2))*(kap^(j))*(factorial(j)/factorial(2))*ET_exp(3,j)*Cp(3,j);
    Dp(3,4)=Dp(3,4)+ ((-1)^(j-3))*(kap^(j))*(factorial(j)/factorial(3))*ET_exp(4,j)*Cp(3,j);
    Dp(3,5)=Dp(3,5)+ ((-1)^(j-4))*(kap^(j))*(factorial(j)/factorial(4))*ET_exp(5,j)*Cp(3,j);
    Dp(3,6)=Dp(3,6)+ ((-1)^(j-5))*(kap^(j))*(factorial(j)/factorial(5))*ET_exp(6,j)*Cp(3,j);
    Dp(3,7)=Dp(3,7)+ ((-1)^(j-6))*(kap^(j))*(factorial(j)/factorial(6))*ET_exp(7,j)*Cp(3,j);
    
end 

% integration constants of main functions
fun_curve_cons_0_exp= s0.*(Dp(1,1)+ Dp(1,2)*log(s0)+ Dp(1,3)*log(s0).^2 + Dp(1,4)*log(s0).^3 + Dp(1,5)*log(s0).^4 + + Dp(1,6)*log(s0).^5 + Dp(1,7)*log(s0).^6);
fun_curve_cons_1_exp= s0.*(Dp(2,1)+ Dp(2,2)*log(s0)+ Dp(2,3)*log(s0).^2 + Dp(2,4)*log(s0).^3 + Dp(2,5)*log(s0).^4 + + Dp(2,6)*log(s0).^5 + Dp(2,7)*log(s0).^6);
fun_curve_cons_2_exp= s0.*(Dp(3,1)+ Dp(3,2)*log(s0)+ Dp(3,3)*(log(s0).^2) + Dp(3,4)*log(s0).^3 + Dp(3,5)*log(s0).^4 + + Dp(3,6)*log(s0).^5 + Dp(3,7)*log(s0).^6);

% main functions 
fun_curve_0_exp= s.*(Dp(1,1)+ Dp(1,2)*log(s)+ Dp(1,3)*log(s).^2 + Dp(1,4)*log(s).^3 + Dp(1,5)*log(s).^4 + + Dp(1,6)*log(s).^5 + Dp(1,7)*log(s).^6) - fun_curve_cons_0_exp;
fun_curve_1_exp= s.*(Dp(2,1)+ Dp(2,2)*log(s)+ Dp(2,3)*log(s).^2 + Dp(2,4)*log(s).^3 + Dp(2,5)*log(s).^4 + + Dp(2,6)*log(s).^5 + Dp(2,7)*log(s).^6) -fun_curve_cons_1_exp;
fun_curve_2_exp= s.*(Dp(3,1)+ Dp(3,2)*log(s)+ Dp(3,3)*(log(s).^2) + Dp(3,4)*log(s).^3 + Dp(3,5)*log(s).^4 + + Dp(3,6)*log(s).^5 + Dp(3,7)*log(s).^6) -fun_curve_cons_2_exp;



% NON-EXPANDED METHOD (not expanding the exponential term in integral to obtain T)

% EN matrix (array with values needed to calculate Cmk values)
EN_nonexp=[[1,0,0,1,1,1];[0,0,0,1,1,1];[0,0,0,1,1,1];[0,0,0,1,1,1];[0,0,0,0,1,1];[0,0,0,0,0,1]];
EN_nonexp(:,:,2)=[[0,1,0,1,1,1];[0,1,0,1,1,1];[0,0,0,1,1,1];[0,0,0,1,1,1];[0,0,0,0,1,1];[0,0,0,0,0,1]];
EN_nonexp(:,:,3)=[[0,0,1,1,1,1];[0,0,1,1,1,1];[0,0,1,1,1,1];[0,0,0,1,1,1];[0,0,0,0,1,1];[0,0,0,0,0,1]];

% Cmk matrix (array of coefficients of F_Tm)
Cmk=zeros(3,6);
for j=0:5
    % Calculating coefficiEN_nonexpts of first function
    Cmk(1,1)=Cmk(1,1) + ((-1)^(j))*factorial(j)*(kap^(j))*EN_nonexp(1,j+1,1)*B(1,j+1);
    Cmk(1,2)=Cmk(1,2) + ((-1)^(j-1))*factorial(j)*(kap^(j-1))*EN_nonexp(2,j+1,1)*B(1,j+1);
    Cmk(1,3)=Cmk(1,3) + ((-1)^(j))*(factorial(j)/factorial(2))*(kap^(j-2))*EN_nonexp(3,j+1,1)*B(1,j+1);
    Cmk(1,4)=Cmk(1,4) + ((-1)^(j-3))*(factorial(j)/factorial(3))*(kap^(j-3))*EN_nonexp(4,j+1,1)*B(1,j+1);
    Cmk(1,5)=Cmk(1,5) + ((-1)^(j-4))*(factorial(j)/factorial(4))*(kap^(j-4))*EN_nonexp(5,j+1,1)*B(1,j+1);
    Cmk(1,6)=Cmk(1,6) + ((-1)^(j-5))*(factorial(j)/factorial(5))*(kap^(j-5))*EN_nonexp(6,j+1,1)*B(1,j+1);
    % Calculating coefficiEN_nonexpts of second function
    Cmk(2,1)=Cmk(2,1) + ((-1)^(j))*factorial(j)*(kap^(j))*EN_nonexp(1,j+1,2)*B(2,j+1);
    Cmk(2,2)=Cmk(2,2) + ((-1)^(j-1))*factorial(j)*(kap^(j-1))*EN_nonexp(2,j+1,2)*B(2,j+1);
    Cmk(2,3)=Cmk(2,3) + ((-1)^(j))*(factorial(j)/factorial(2))*(kap^(j-2))*EN_nonexp(3,j+1,2)*B(2,j+1);
    Cmk(2,4)=Cmk(2,4) + ((-1)^(j-3))*(factorial(j)/factorial(3))*(kap^(j-3))*EN_nonexp(4,j+1,2)*B(2,j+1);
    Cmk(2,5)=Cmk(2,5) + ((-1)^(j-4))*(factorial(j)/factorial(4))*(kap^(j-4))*EN_nonexp(5,j+1,2)*B(2,j+1);
    Cmk(2,6)=Cmk(2,6) + ((-1)^(j-5))*(factorial(j)/factorial(5))*(kap^(j-5))*EN_nonexp(6,j+1,2)*B(2,j+1);
    % Calculating coefficiEN_nonexpts of third function
    Cmk(3,1)=Cmk(3,1) + ((-1)^(j))*factorial(j)*(kap^(j))*EN_nonexp(1,j+1,3)*B(3,j+1);
    Cmk(3,2)=Cmk(3,2) + ((-1)^(j-1))*factorial(j)*(kap^(j-1))*EN_nonexp(2,j+1,3)*B(3,j+1);
    Cmk(3,3)=Cmk(3,3) + ((-1)^(j))*(factorial(j)/factorial(2))*(kap^(j-2))*EN_nonexp(3,j+1,3)*B(3,j+1);
    Cmk(3,4)=Cmk(3,4) + ((-1)^(j-3))*(factorial(j)/factorial(3))*(kap^(j-3))*EN_nonexp(4,j+1,3)*B(3,j+1);
    Cmk(3,5)=Cmk(3,5) + ((-1)^(j-4))*(factorial(j)/factorial(4))*(kap^(j-4))*EN_nonexp(5,j+1,3)*B(3,j+1);
    Cmk(3,6)=Cmk(3,6) + ((-1)^(j-5))*(factorial(j)/factorial(5))*(kap^(j-5))*EN_nonexp(6,j+1,3)*B(3,j+1);
end

% ET matrix (array with values needed to calculate Dqk values)
ET_nonexp=[[1/2, 1/4, 1/4, 3/8, 3/4, 15/8];[0, 1/2, 1/2, 3/4, 3/2, 15/4];[0, 0, 1/2, 3/4, 3/2, 15/4];[0, 0, 0, 1/2, 1, 5/2];[0, 0, 0, 0, 1/2, 5/4];[0, 0, 0, 0, 0, 1/2]];

% Dqk matrix (array of coefficients of curve function terms)
Dqk=zeros(3,6);
for j=0:5
    % Calculating coefficients of first function
    Dqk(1,1)=Dqk(1,1)+ ((-1)^(j))*(kap^(j))*ET_nonexp(1,j+1)*Cmk(1,j+1);
    Dqk(1,2)=Dqk(1,2)+ ((-1)^(j-1))*(kap^(j))*ET_nonexp(2,j+1)*Cmk(1,j+1);
    Dqk(1,3)=Dqk(1,3)+ ((-1)^(j-2))*(kap^(j))*ET_nonexp(3,j+1)*Cmk(1,j+1);
    Dqk(1,4)=Dqk(1,4)+ ((-1)^(j-3))*(kap^(j))*ET_nonexp(4,j+1)*Cmk(1,j+1);
    Dqk(1,5)=Dqk(1,5)+ ((-1)^(j-4))*(kap^(j))*ET_nonexp(5,j+1)*Cmk(1,j+1);
    Dqk(1,6)=Dqk(1,6)+ ((-1)^(j-5))*(kap^(j))*ET_nonexp(6,j+1)*Cmk(1,j+1);
    % Calculating coefficients of second function
    Dqk(2,1)=Dqk(2,1)+ ((-1)^(j))*(kap^(j))*ET_nonexp(1,j+1)*Cmk(2,j+1);
    Dqk(2,2)=Dqk(2,2)+ ((-1)^(j-1))*(kap^(j))*ET_nonexp(2,j+1)*Cmk(2,j+1);
    Dqk(2,3)=Dqk(2,3)+ ((-1)^(j-2))*(kap^(j))*ET_nonexp(3,j+1)*Cmk(2,j+1);
    Dqk(2,4)=Dqk(2,4)+ ((-1)^(j-3))*(kap^(j))*ET_nonexp(4,j+1)*Cmk(2,j+1);
    Dqk(2,5)=Dqk(2,5)+ ((-1)^(j-4))*(kap^(j))*ET_nonexp(5,j+1)*Cmk(2,j+1);
    Dqk(2,6)=Dqk(2,6)+ ((-1)^(j-5))*(kap^(j))*ET_nonexp(6,j+1)*Cmk(2,j+1);
    % Calculating coefficients of third function
    Dqk(3,1)=Dqk(3,1)+ ((-1)^(j))*(kap^(j))*ET_nonexp(1,j+1)*Cmk(3,j+1);
    Dqk(3,2)=Dqk(3,2)+ ((-1)^(j-1))*(kap^(j))*ET_nonexp(2,j+1)*Cmk(3,j+1);
    Dqk(3,3)=Dqk(3,3)+ ((-1)^(j-2))*(kap^(j))*ET_nonexp(3,j+1)*Cmk(3,j+1);
    Dqk(3,4)=Dqk(3,4)+ ((-1)^(j-3))*(kap^(j))*ET_nonexp(4,j+1)*Cmk(3,j+1);
    Dqk(3,5)=Dqk(3,5)+ ((-1)^(j-4))*(kap^(j))*ET_nonexp(5,j+1)*Cmk(3,j+1);
    Dqk(3,6)=Dqk(3,6)+ ((-1)^(j-5))*(kap^(j))*ET_nonexp(6,j+1)*Cmk(3,j+1);
end



% CONSTRUCTING FUNCTIONS, VECTORS, AND CURVE AFTER CALCULATING COEFFICIENTS

% integration constants of main functions
fun_curve_cons_0_nonexp=kap*(s0.^2).*(Dqk(1,1)+ Dqk(1,2)*log(s0)+ Dqk(1,3)*(log(s0).^2) + Dqk(1,4)*(log(s0).^3) + Dqk(1,5)*(log(s0).^4) + + Dqk(1,6)*(log(s0).^5));
fun_curve_cons_1_nonexp=kap*(s0.^2).*(Dqk(2,1)+ Dqk(2,2)*log(s0)+ Dqk(2,3)*(log(s0).^2) + Dqk(2,4)*(log(s0).^3) + Dqk(2,5)*(log(s0).^4) + + Dqk(2,6)*(log(s0).^5));
fun_curve_cons_2_nonexp=kap*(s0.^2).*(Dqk(3,1)+ Dqk(3,2)*log(s0)+ Dqk(3,3)*(log(s0).^2) + Dqk(3,4)*(log(s0).^3) + Dqk(3,5)*(log(s0).^4) + + Dqk(3,6)*(log(s0).^5));

% main functions
fun_curve_0_nonexp= kap*(s.^2).*(Dqk(1,1)+ Dqk(1,2)*log(s)+ Dqk(1,3)*(log(s).^2) + Dqk(1,4)*(log(s).^3) + Dqk(1,5)*(log(s).^4) + + Dqk(1,6)*(log(s).^5)) - kap*Cmk(1,1)*(s-s0)  - fun_curve_cons_0_nonexp;
fun_curve_1_nonexp= kap*(s.^2).*(Dqk(2,1)+ Dqk(2,2)*log(s)+ Dqk(2,3)*(log(s).^2) + Dqk(2,4)*(log(s).^3) + Dqk(2,5)*(log(s).^4) + + Dqk(2,6)*(log(s).^5)) - kap*Cmk(2,1)*(s-s0) - fun_curve_cons_1_nonexp;
fun_curve_2_nonexp= kap*(s.^2).*(Dqk(3,1)+ Dqk(3,2)*log(s)+ Dqk(3,3)*((log(s).^2)) + Dqk(3,4)*(log(s).^3) + Dqk(3,5)*(log(s).^4) + + Dqk(3,6)*(log(s).^5)) - kap*Cmk(3,1)*(s-s0) - fun_curve_cons_2_nonexp;


% FINAL EQUATION OF CURVES (EXPANDED AND NON-EXPANDED METHODS)
curve_exp(1,:)= s-s0-(1/(2*kap))*fun_curve_2_exp - fun_curve_1_exp; %COEFFICIENT OF VECTOR A
curve_exp(2,:)=fun_curve_0_exp-fun_curve_2_exp; %COEFFICIENT OF VECTOR B
curve_exp(3,:)=-fun_curve_1_exp; %COEFFICIENT OF VECTOR C
curve_nonexp(1,:)= s-s0-(1/(2*kap))*fun_curve_2_nonexp - fun_curve_1_nonexp; %COEFFICIENT OF VECTOR A
curve_nonexp(2,:)=fun_curve_0_nonexp-fun_curve_2_nonexp;    %COEFFICIENT OF VECTOR B
curve_nonexp(3,:)=-fun_curve_1_nonexp; %COEFFICIENT OF VECTOR C


% functions needed for normal vector
fun_norm_0= B(1,1) + B(1,2)*t + + B(1,3)*t.^2 + B(1,4)*t.^3 + B(1,5)*t.^4 + B(1,6)*t.^5;
fun_norm_1= B(2,1) + B(2,2)*t + + B(2,3)*t.^2 + B(2,4)*t.^3 + B(2,5)*t.^4 + B(2,6)*t.^5;
fun_norm_2= B(3,1) + B(3,2)*t + + B(3,3)*t.^2 + B(3,4)*t.^3 + B(3,5)*t.^4 + B(3,6)*t.^5;

% PRINCIPAL NORMAL VECTOR (SAME FOR BOTH METHODS)
norm(1,:)=-fun_norm_1 - (1/(2*kap))*fun_norm_2;
norm(2,:)= fun_norm_0 - fun_norm_2;
norm(3,:)=fun_norm_1;

% functions needed for tangent vector (non-expanded method)
fun_tang_nonexp_0= kap*exp(t./kap).*(Cmk(1,1) + Cmk(1,2)*t + + Cmk(1,3)*t.^2 + Cmk(1,4)*t.^3 + Cmk(1,5)*t.^4 + Cmk(1,6)*t.^5) - kap*Cmk(1,1);
fun_tang_nonexp_1= kap*exp(t./kap).*(Cmk(2,1) + Cmk(2,2)*t + + Cmk(2,3)*t.^2 + Cmk(2,4)*t.^3 + Cmk(2,5)*t.^4 + Cmk(2,6)*t.^5) - kap*Cmk(2,1);
fun_tang_nonexp_2= kap*exp(t./kap).*(Cmk(3,1) + Cmk(3,2)*t + + Cmk(3,3)*t.^2 + Cmk(3,4)*t.^3 + Cmk(3,5)*t.^4 + Cmk(3,6)*t.^5) - kap*Cmk(3,1);

% TANGENT VECTOR (NON-EXPANDED METHOD)
tang_nonexp(1,:)=1-fun_tang_nonexp_1 - (1/(2*kap))*fun_tang_nonexp_2;
tang_nonexp(2,:)= fun_tang_nonexp_0 - fun_tang_nonexp_2;
tang_nonexp(3,:)=fun_tang_nonexp_1;

% BINORMAL VECTOR (NON-EXPANDED METHOD)
bino_nonexp(1,:)=tang_nonexp(2,:).*norm(3,:)-tang_nonexp(3,:).*norm(2,:);
bino_nonexp(2,:)=tang_nonexp(3,:).*norm(1,:)-tang_nonexp(1,:).*norm(3,:);
bino_nonexp(3,:)=tang_nonexp(1,:).*norm(2,:)-tang_nonexp(2,:).*norm(1,:);

% functions needed for tangent vector
fun_tang_exp_0= Cp(1,1)*t + Cp(1,2)*t.^2 + + Cp(1,3)*t.^3 + Cp(1,4)*t.^4 + Cp(1,5)*t.^5 + Cp(1,6)*t.^6;
fun_tang_exp_1= Cp(2,1)*t + Cp(2,2)*t.^2 + + Cp(2,3)*t.^3 + Cp(2,4)*t.^4 + Cp(2,5)*t.^5 + Cp(2,6)*t.^6;
fun_tang_exp_2= Cp(3,1)*t + Cp(3,2)*t.^2 + + Cp(3,3)*t.^3 + Cp(3,4)*t.^4 + Cp(3,5)*t.^5 + Cp(3,6)*t.^6;

% TANGENT VECTOR (EXPANDED METHOD)
tang_exp(1,:)=1-fun_tang_exp_1 - (1/(2*kap))*fun_tang_exp_2;
tang_exp(2,:)= fun_tang_exp_0 - fun_tang_exp_2;
tang_exp(3,:)=fun_tang_exp_1;

% CALCULATING MAGNITUDE AND ORTHOGONALITY OF FRENET FRAME VECTORS

% BINORMAL VECTOR (NON-EXPANDED METHOD)
bino_exp(1,:)=tang_exp(2,:).*norm(3,:)-tang_exp(3,:).*norm(2,:);
bino_exp(2,:)=tang_exp(3,:).*norm(1,:)-tang_exp(1,:).*norm(3,:);
bino_exp(3,:)=tang_exp(1,:).*norm(2,:)-tang_exp(2,:).*norm(1,:);

% "s" CALCULATING NORMS OF FRENET FRAME FOR EACH "s"
magn(1,:)=tang_exp(1,:).^2 + tang_exp(2,:).^3 + tang_exp(2,:).^3;
magn(2,:)=tang_nonexp(1,:).^2 + tang_nonexp(2,:).^3 + tang_nonexp(2,:).^3;
magn(3,:)=norm(1,:).^2 + norm(2,:).^3 + norm(3,:).^3;
magn(4,:)=bino_exp(1,:).^2 + bino_exp(2,:).^3 + bino_exp(3,:).^3;
magn(5,:)=bino_nonexp(1,:).^2 + bino_nonexp(2,:).^3 + bino_nonexp(3,:).^3;

% DOT PRODUCTS BETWEEN FRENET FRAME VECTORS FOR EACH "s"
dot_tnb(1,:)=tang_nonexp(1,:).*norm(1,:) + tang_nonexp(2,:).*norm(2,:) + tang_nonexp(3,:).*norm(3,:);
dot_tnb(2,:)=tang_nonexp(1,:).*bino_nonexp(1,:) + tang_nonexp(2,:).*bino_nonexp(2,:) + tang_nonexp(3,:).*bino_nonexp(3,:);
dot_tnb(3,:)=norm(1,:).*bino_nonexp(1,:) + norm(2,:).*bino_nonexp(2,:) + norm(3,:).*bino_nonexp(3,:);
dot_tnb(4,:)=tang_exp(1,:).*norm(1,:) + tang_exp(2,:).*norm(2,:) + tang_exp(3,:).*norm(3,:);
dot_tnb(5,:)=tang_exp(1,:).*bino_exp(1,:) + tang_exp(2,:).*bino_exp(2,:) + tang_exp(3,:).*bino_exp(3,:);
dot_tnb(6,:)=norm(1,:).*bino_exp(1,:) + norm(2,:).*bino_exp(2,:) + norm(3,:).*bino_exp(3,:);

% EXTRA DATA
% derivative of functions for normal vector
dfunnorm_dt_0=3*B(1,4)*t.^2 + 4*B(1,5)*t.^3 + 5*B(1,6)*t.^4;
dfunnorm_dt_1=B(2,2)+3*B(2,4)*t.^2 + 4*B(2,5)*t.^3 + 5*B(2,6)*t.^4;
dfunnorm_dt_2=2*B(3,3)*t+3*B(3,4)*t.^2 + 4*B(3,5)*t.^3 + 5*B(3,6)*t.^4;

% derivative of normal vector
dnorm_dt(1,:)=-dfunnorm_dt_1 - (1/(2*kap))*dfunnorm_dt_2;
dnorm_dt(2,:)= dfunnorm_dt_0 - dfunnorm_dt_2;
dnorm_dt(3,:)=dfunnorm_dt_1;

% binormal vector from Eq. 3.9
bino_exp_2(1,:)=dnorm_dt(1,:)+exp(t./kap).*tang_exp(1,:);
bino_exp_2(2,:)=dnorm_dt(2,:)+exp(t./kap).*tang_exp(2,:);
bino_exp_2(3,:)=dnorm_dt(3,:)+exp(t./kap).*tang_exp(3,:);
bino_nonexp_2(1,:)=dnorm_dt(1,:)+exp(t./kap).*tang_nonexp(1,:);
bino_nonexp_2(2,:)=dnorm_dt(2,:)+exp(t./kap).*tang_nonexp(2,:);
bino_nonexp_2(3,:)=dnorm_dt(3,:)+exp(t./kap).*tang_nonexp(3,:);

% magnitude of binormal vector
magn_bino(1,:)=bino_exp(1,:).^2 + bino_exp(2,:).^3 + bino_exp(3,:).^3;
magn_bino(2,:)=bino_exp_2(1,:).^2 + bino_exp_2(2,:).^3 + bino_exp_2(3,:).^3;
magn_bino(3,:)=bino_nonexp(1,:).^2 + bino_nonexp(2,:).^3 + bino_nonexp(3,:).^3;
magn_bino(4,:)=bino_nonexp_2(1,:).^2 + bino_nonexp_2(2,:).^3 + bino_nonexp_2(3,:).^3;

figure(1)
subplot(1,2,1)
plot(s,magn_bino(1,:))
hold on
plot(s,magn_bino(2,:))
subplot(1,2,2)
plot(s,magn_bino(3,:))
hold on
plot(s,magn_bino(4,:))


%{
% PLOTTING DATA

% PLOTTING CURVE (EXPANDED AND NON-EXPANDED)
figure(1)
camproj('orthographic')
hold off
plot3(curve_exp(1,:), curve_exp(2,:), curve_exp(3,:), "b")
axis equal
grid on
hold on
plot3(curve_nonexp(1,:), curve_nonexp(2,:), curve_nonexp(3,:), "r")
legend("Expanded method","Non-Expanded method")
xlabel("x")
ylabel("y")
zlabel("z")

% PLOTTING F_cm f FOR EXPANDED AND NON-EXPANDED METHODS
figure(2)
subplot(1,3,1)
plot(s,fun_curve_0_exp, "r") % f_c0 function obtained from expanded method
hold on
plot(s,fun_curve_0_nonexp, "b") % f_c0 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C0}(s)")
legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(1,3,2)
plot(s,fun_curve_1_exp, "r") % f_c1 function obtained from expanded method
hold on
plot(s,fun_curve_1_nonexp, "b") % f_c1 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C1}(s)")
legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(1,3,3)
plot(s,fun_curve_2_exp, "r") % f_c2 function obtained from expanded method
hold on
plot(s,fun_curve_2_nonexp, "b") % f_c2 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C2}(s)")
legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off


% PLOTTING MAGNITUDE OF VECTORS AND DOT PRODUCTS FOR EXPANDED AND
% NON-EXPANDED METHODS
figure(3)
subplot(2,2,1)
plot(s,magn(1,:),'r')
hold on
plot(s,magn(3,:),'b')
plot(s,magn(4,:),'k')
xlabel("s")
ylabel("|X|")
title("Expanded Method (\kappa=" + kap + ")")
legend("Tangent Vector","Principal Normal Vector", "Binormal Vector","Location","Northwest")
hold off
subplot(2,2,2)
plot(s,dot_tnb(4,:),'r')
hold on
plot(s,dot_tnb(5,:),'b')
plot(s,dot_tnb(6,:),'k')
xlabel("s")
ylabel("X_{1}\cdotX_{2}")
title("Expanded Method (\kappa=" + kap + ")")
legend("T \cdot N", "T \cdot B", "N \cdot B","Location","Southwest")
hold off
subplot(2,2,3)
plot(s,magn(2,:),'r')
hold on
plot(s,magn(3,:),'b')
plot(s,magn(5,:),'k')
xlabel("s")
ylabel("|X|")
title("Non-Expanded Method (\kappa=" + kap + ")")
legend("Tangent Vector","Principal Normal Vector", "Binormal Vector","Location","Southwest")
hold off
subplot(2,2,4)
plot(s,dot_tnb(1,:),'r')
hold on
plot(s,dot_tnb(2,:),'b')
plot(s,dot_tnb(3,:),'k')
xlabel("s")
ylabel("X_{1}\cdotX_{2}")
title("Non-Expanded Method (\kappa=" + kap + ")")
legend("T \cdot N", "T \cdot B", "N \cdot B","Location","Southwest")
hold off


% PLOTTING F_tm AND F_cm f FOR EXPANDED AND NON-EXPANDED METHODS
figure(4)
subplot(2,3,1)
plot(s,fun_tang_exp_0, "r") % f_t0 function obtained from expanded method
hold on
plot(s,fun_tang_nonexp_0, "b") % f_t0 function obtained from non-expanded method
xlabel("s")
ylabel("f_{T0}(s)")
legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(2,3,2)
plot(s,fun_tang_exp_1, "r") % f_t1 function obtained from expanded method
hold on
plot(s,fun_tang_nonexp_1, "b") % f_t1 function obtained from non-expanded method
xlabel("s")
ylabel("f_{T1}(s)")
legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(2,3,3)
plot(s,fun_tang_exp_0, "r") % f_t2 function obtained from expanded method
hold on
plot(s,fun_tang_nonexp_0, "b") % f_t2 function obtained from non-expanded method
xlabel("s")
ylabel("f_{T2}(s)")
legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(2,3,4)
plot(s,fun_curve_0_exp, "r") % f_c0 function obtained from expanded method
hold on
plot(s,fun_curve_0_nonexp, "b") % f_c0 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C0}(s)")
legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(2,3,5)
plot(s,fun_curve_1_exp, "r") % f_c1 function obtained from expanded method
hold on
plot(s,fun_curve_1_nonexp, "b") % f_c1 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C1}(s)")
legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(2,3,6)
plot(s,fun_curve_2_exp, "r") % f_c2 function obtained from expanded method
hold on
plot(s,fun_curve_2_nonexp, "b") % f_c2 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C2}(s)")
legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off

%}























































































































