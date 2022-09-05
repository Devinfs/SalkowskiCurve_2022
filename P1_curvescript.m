% Code needed to caluclate and plot quantities
% from expanded and non-expanded methods.
% Quantities that should be changed will have (CHANGE) next to them.

clear


% GENERAL PARAMETERS/INFORMATION
% arc length domain
s0=1; % initial arclength (do not change)
sf=5; % final arclength value (CHANGE)
ds=0.00002; % step size (CHANGE)
s=s0:ds:sf;
%s=s0:0.000000002:1.0001; used for showing how small interval around s=1


% value of curvature
kap=1; %(CHANGE)

% t interval
t=kap*log(s);



% tangent vector coefficient functions (F_Tm)
% normal vector coefficient functions (F_Nm)

% B matrix (array of coefficients of F_Nm)(EQ. 4.5)
B=[[1,0,0, -1/(factorial(3)*kap),-1/(factorial(4)*(kap^2)),(1/factorial(5))*((-13/(kap^3))+(1/kap))];[0,1,0, -2/factorial(3),-5/(factorial(4)*kap),(1/factorial(5))*((-17/(kap^2))+4)];[0,0,1,2/(factorial(3)*kap),(1/factorial(4))*((2/(kap^2))-4),(1/factorial(5))*((2/(kap^3))-(18/kap))]];


% EXPANDED METHOD (expanding the exponential term in integral to obtain T)

% G_mp matrix (array of coeffecients of F_Tm) (EQ. 5.4)
G_mp= [[1, 1/(2*kap), 1/(factorial(3)*kap^2), (1/(factorial(4)))*((1/kap^3)-(1/kap)), (1/factorial(5))*((1/kap^4)-(9/kap^2)),0,0,0,0,0]; [0, 1/2, 1/(3*kap), ((1/(8*kap^2))-(1/12)),((1/(30*kap^3))-(13/(120*kap^2))),0,0,0,0,0];[0,0,1/3,1/(3*kap),((-1/30)+(11/(60*kap^2))), ((1/(45*kap^3))-(19/(360*kap)))],0,0,0,0];


% The next two "%" relate to calculating G_mp coefficients from Appendix A
% ET matrix (array with values needed to calculate G_mp values)
ET_exp=[[1, 1, 1, 1, 1, 1]; [1, 1, 1, 1, 1, 1]; [0, 1, 1, 1, 1, 1]; [0, 0, 1, 1, 1, 1]; [0, 0, 0, 1, 1, 1]; [0, 0, 0, 0, 1, 1]; [0, 0, 0, 0, 0, 1]];

% H_mq matrix (array of coefficients of main functions)
H_mq=zeros(3,7);
for j=1:6
    % Calculating coefficients of first function
    H_mq(1,1)=H_mq(1,1)+ ((-1)^(j))*(kap^(j))*factorial(j)*ET_exp(1,j)*G_mp(1,j);
    H_mq(1,2)=H_mq(1,2)+ ((-1)^(j-1))*(kap^(j))*factorial(j)*ET_exp(2,j)*G_mp(1,j);
    H_mq(1,3)=H_mq(1,3)+ ((-1)^(j-2))*(kap^(j))*(factorial(j)/factorial(2))*ET_exp(3,j)*G_mp(1,j);
    H_mq(1,4)=H_mq(1,4)+ ((-1)^(j-3))*(kap^(j))*(factorial(j)/factorial(3))*ET_exp(4,j)*G_mp(1,j);
    H_mq(1,5)=H_mq(1,5)+ ((-1)^(j-4))*(kap^(j))*(factorial(j)/factorial(4))*ET_exp(5,j)*G_mp(1,j);
    H_mq(1,6)=H_mq(1,6)+ ((-1)^(j-5))*(kap^(j))*(factorial(j)/factorial(5))*ET_exp(6,j)*G_mp(1,j);
    H_mq(1,7)=H_mq(1,7)+ ((-1)^(j-6))*(kap^(j))*(factorial(j)/factorial(6))*ET_exp(7,j)*G_mp(1,j);
    % Calculating coefficients of second function
    H_mq(2,1)=H_mq(2,1)+ ((-1)^(j))*(kap^(j))*factorial(j)*ET_exp(1,j)*G_mp(2,j);
    H_mq(2,2)=H_mq(2,2)+ ((-1)^(j-1))*(kap^(j))*factorial(j)*ET_exp(2,j)*G_mp(2,j);
    H_mq(2,3)=H_mq(2,3)+ ((-1)^(j-2))*(kap^(j))*(factorial(j)/factorial(2))*ET_exp(3,j)*G_mp(2,j);
    H_mq(2,4)=H_mq(2,4)+ ((-1)^(j-3))*(kap^(j))*(factorial(j)/factorial(3))*ET_exp(4,j)*G_mp(2,j);
    H_mq(2,5)=H_mq(2,5)+ ((-1)^(j-4))*(kap^(j))*(factorial(j)/factorial(4))*ET_exp(5,j)*G_mp(2,j);
    H_mq(2,6)=H_mq(2,6)+ ((-1)^(j-5))*(kap^(j))*(factorial(j)/factorial(5))*ET_exp(6,j)*G_mp(2,j);
    H_mq(2,7)=H_mq(2,7)+ ((-1)^(j-6))*(kap^(j))*(factorial(j)/factorial(6))*ET_exp(7,j)*G_mp(2,j);
    % Calculating coefficients of third function
    H_mq(3,1)=H_mq(3,1)+ ((-1)^(j))*(kap^(j))*factorial(j)*ET_exp(1,j)*G_mp(3,j);
    H_mq(3,2)=H_mq(3,2)+ ((-1)^(j-1))*(kap^(j))*factorial(j)*ET_exp(2,j)*G_mp(3,j);
    H_mq(3,3)=H_mq(3,3)+ ((-1)^(j-2))*(kap^(j))*(factorial(j)/factorial(2))*ET_exp(3,j)*G_mp(3,j);
    H_mq(3,4)=H_mq(3,4)+ ((-1)^(j-3))*(kap^(j))*(factorial(j)/factorial(3))*ET_exp(4,j)*G_mp(3,j);
    H_mq(3,5)=H_mq(3,5)+ ((-1)^(j-4))*(kap^(j))*(factorial(j)/factorial(4))*ET_exp(5,j)*G_mp(3,j);
    H_mq(3,6)=H_mq(3,6)+ ((-1)^(j-5))*(kap^(j))*(factorial(j)/factorial(5))*ET_exp(6,j)*G_mp(3,j);
    H_mq(3,7)=H_mq(3,7)+ ((-1)^(j-6))*(kap^(j))*(factorial(j)/factorial(6))*ET_exp(7,j)*G_mp(3,j);
    
end 

% The next two "%" relate to the final functions needed for curve (EQ. 8.1)
% integration constants of main functions
fun_curve_cons_0_exp= s0.*(H_mq(1,1)+ H_mq(1,2)*log(s0)+ H_mq(1,3)*log(s0).^2 + H_mq(1,4)*log(s0).^3 + H_mq(1,5)*log(s0).^4 + + H_mq(1,6)*log(s0).^5 + H_mq(1,7)*log(s0).^6);
fun_curve_cons_1_exp= s0.*(H_mq(2,1)+ H_mq(2,2)*log(s0)+ H_mq(2,3)*log(s0).^2 + H_mq(2,4)*log(s0).^3 + H_mq(2,5)*log(s0).^4 + + H_mq(2,6)*log(s0).^5 + H_mq(2,7)*log(s0).^6);
fun_curve_cons_2_exp= s0.*(H_mq(3,1)+ H_mq(3,2)*log(s0)+ H_mq(3,3)*(log(s0).^2) + H_mq(3,4)*log(s0).^3 + H_mq(3,5)*log(s0).^4 + + H_mq(3,6)*log(s0).^5 + H_mq(3,7)*log(s0).^6);

% main functions 
fun_curve_0_exp= s.*(H_mq(1,1)+ H_mq(1,2)*log(s)+ H_mq(1,3)*log(s).^2 + H_mq(1,4)*log(s).^3 + H_mq(1,5)*log(s).^4 + + H_mq(1,6)*log(s).^5 + H_mq(1,7)*log(s).^6) - fun_curve_cons_0_exp;
fun_curve_1_exp= s.*(H_mq(2,1)+ H_mq(2,2)*log(s)+ H_mq(2,3)*log(s).^2 + H_mq(2,4)*log(s).^3 + H_mq(2,5)*log(s).^4 + + H_mq(2,6)*log(s).^5 + H_mq(2,7)*log(s).^6) -fun_curve_cons_1_exp;
fun_curve_2_exp= s.*(H_mq(3,1)+ H_mq(3,2)*log(s)+ H_mq(3,3)*(log(s).^2) + H_mq(3,4)*log(s).^3 + H_mq(3,5)*log(s).^4 + + H_mq(3,6)*log(s).^5 + H_mq(3,7)*log(s).^6) -fun_curve_cons_2_exp;



% NON-EXPANDED METHOD (not expanding the exponential term in integral to obtain T)

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
fun_curve_cons_0_nonexp=kap*(s0.^2).*(D_ml(1,1)+ D_ml(1,2)*log(s0)+ D_ml(1,3)*(log(s0).^2) + D_ml(1,4)*(log(s0).^3) + D_ml(1,5)*(log(s0).^4) + + D_ml(1,6)*(log(s0).^5));
fun_curve_cons_1_nonexp=kap*(s0.^2).*(D_ml(2,1)+ D_ml(2,2)*log(s0)+ D_ml(2,3)*(log(s0).^2) + D_ml(2,4)*(log(s0).^3) + D_ml(2,5)*(log(s0).^4) + + D_ml(2,6)*(log(s0).^5));
fun_curve_cons_2_nonexp=kap*(s0.^2).*(D_ml(3,1)+ D_ml(3,2)*log(s0)+ D_ml(3,3)*(log(s0).^2) + D_ml(3,4)*(log(s0).^3) + D_ml(3,5)*(log(s0).^4) + + D_ml(3,6)*(log(s0).^5));

% main functions
fun_curve_0_nonexp= kap*(s.^2).*(D_ml(1,1)+ D_ml(1,2)*log(s)+ D_ml(1,3)*(log(s).^2) + D_ml(1,4)*(log(s).^3) + D_ml(1,5)*(log(s).^4) + + D_ml(1,6)*(log(s).^5)) - fun_curve_cons_0_nonexp - kap*Cmk(1,1)*(s-s0);
fun_curve_1_nonexp= kap*(s.^2).*(D_ml(2,1)+ D_ml(2,2)*log(s)+ D_ml(2,3)*(log(s).^2) + D_ml(2,4)*(log(s).^3) + D_ml(2,5)*(log(s).^4) + + D_ml(2,6)*(log(s).^5)) - fun_curve_cons_1_nonexp - kap*Cmk(2,1)*(s-s0);
fun_curve_2_nonexp= kap*(s.^2).*(D_ml(3,1)+ D_ml(3,2)*log(s)+ D_ml(3,3)*((log(s).^2)) + D_ml(3,4)*(log(s).^3) + D_ml(3,5)*(log(s).^4) + + D_ml(3,6)*(log(s).^5)) - fun_curve_cons_2_nonexp - kap*Cmk(3,1)*(s-s0);


% FINAL EQUATION OF CURVES (EXPANDED AND NON-EXPANDED METHODS)(THEOREM 3&4)
% EXPANDED METHOD
curve_exp(1,:)= s-s0-(1/(2*kap))*fun_curve_2_exp - fun_curve_1_exp; %COEFFICIENT OF VECTOR A
curve_exp(2,:)=fun_curve_0_exp-fun_curve_2_exp; %COEFFICIENT OF VECTOR B
curve_exp(3,:)=-fun_curve_1_exp; %COEFFICIENT OF VECTOR C
% NON-EXPANDED METHOD
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

% functions needed for tangent vector (expanded)
fun_tang_exp_0= G_mp(1,1)*t + G_mp(1,2)*t.^2 + + G_mp(1,3)*t.^3 + G_mp(1,4)*t.^4 + G_mp(1,5)*t.^5 + G_mp(1,6)*t.^6;
fun_tang_exp_1= G_mp(2,1)*t + G_mp(2,2)*t.^2 + + G_mp(2,3)*t.^3 + G_mp(2,4)*t.^4 + G_mp(2,5)*t.^5 + G_mp(2,6)*t.^6;
fun_tang_exp_2= G_mp(3,1)*t + G_mp(3,2)*t.^2 + + G_mp(3,3)*t.^3 + G_mp(3,4)*t.^4 + G_mp(3,5)*t.^5 + G_mp(3,6)*t.^6;

% TANGENT VECTOR (EXPANDED METHOD)
tang_exp(1,:)=1-fun_tang_exp_1 - (1/(2*kap))*fun_tang_exp_2;
tang_exp(2,:)= fun_tang_exp_0 - fun_tang_exp_2;
tang_exp(3,:)=fun_tang_exp_1;

% BINORMAL VECTOR (NON-EXPANDED METHOD)
bino_exp(1,:)=tang_exp(2,:).*norm(3,:)-tang_exp(3,:).*norm(2,:);
bino_exp(2,:)=tang_exp(3,:).*norm(1,:)-tang_exp(1,:).*norm(3,:);
bino_exp(3,:)=tang_exp(1,:).*norm(2,:)-tang_exp(2,:).*norm(1,:);


% CALCULATING MAGNITUDE AND ORTHOGONALITY OF FRENET FRAME VECTORS

% CALCULATING NORMS OF FRENET FRAME FOR EACH "s"
magn(1,:)=sqrt(tang_exp(1,:).^2 + tang_exp(2,:).^2 + tang_exp(3,:).^2);
magn(2,:)=sqrt(tang_nonexp(1,:).^2 + tang_nonexp(2,:).^2 + tang_nonexp(3,:).^2);
magn(3,:)=sqrt(norm(1,:).^2 + norm(2,:).^2 + norm(3,:).^2);
magn(4,:)=sqrt(bino_exp(1,:).^2 + bino_exp(2,:).^2 + bino_exp(3,:).^2);
magn(5,:)=sqrt(bino_nonexp(1,:).^2 + bino_nonexp(2,:).^2 + bino_nonexp(3,:).^2);

% DOT PRODUCTS BETWEEN FRENET FRAME VECTORS FOR EACH "s"
dot_tnb(1,:)=tang_exp(1,:).*norm(1,:) + tang_exp(2,:).*norm(2,:) + tang_exp(3,:).*norm(3,:);
dot_tnb(2,:)=tang_exp(1,:).*bino_exp(1,:) + tang_exp(2,:).*bino_exp(2,:) + tang_exp(3,:).*bino_exp(3,:);
dot_tnb(3,:)=norm(1,:).*bino_exp(1,:) + norm(2,:).*bino_exp(2,:) + norm(3,:).*bino_exp(3,:);
dot_tnb(4,:)=tang_nonexp(1,:).*norm(1,:) + tang_nonexp(2,:).*norm(2,:) + tang_nonexp(3,:).*norm(3,:);
dot_tnb(5,:)=tang_nonexp(1,:).*bino_nonexp(1,:) + tang_nonexp(2,:).*bino_nonexp(2,:) + tang_nonexp(3,:).*bino_nonexp(3,:);
dot_tnb(6,:)=norm(1,:).*bino_nonexp(1,:) + norm(2,:).*bino_nonexp(2,:) + norm(3,:).*bino_nonexp(3,:);


% PLOTTING DATA


%
% FIGURE 4 IN PAPER
% PLOTTING CURVE (EXPANDED AND NON-EXPANDED)
figure(1)
camproj('orthographic')
hold off
plot3(curve_exp(1,:), curve_exp(2,:), curve_exp(3,:), "r",'linewidth',2)
axis equal
grid on
hold on
plot3(curve_nonexp(1,:), curve_nonexp(2,:), curve_nonexp(3,:), "b",'linewidth',2)
plot3(curve_exp(1,1), curve_exp(2,1), curve_exp(3,1), '-o','Color','b','MarkerSize',8,'MarkerFaceColor','#D9FFFF')
xlabel("x")
ylabel("y")
zlabel("z")


%{

% FIGURES 2&3 IN PAPER
% PLOTTING CURVE (EXPANDED/NON-EXPANDED CURVE SOLO)
figure(2)
camproj('orthographic')
hold off
plot3(curve_exp(1,:), curve_exp(2,:), curve_exp(3,:), "B",'linewidth',2)
axis equal
grid on
hold on
plot3(curve_exp(1,1), curve_exp(2,1), curve_exp(3,1), '-o','Color','b','MarkerSize',8,'MarkerFaceColor','#D9FFFF')
xlabel("x")
ylabel("y")
zlabel("z")

%}


%{
% FIGURE 5 IN PAPER
% Calculating difference in position 
diff(1,:)=curve_nonexp(1,:)-curve_exp(1,:);
diff(2,:)=curve_nonexp(2,:)-curve_exp(2,:);
diff(3,:)=curve_nonexp(3,:)-curve_exp(3,:);

figure(3)
plot(s,diff(1,:),"r",'linewidth',2)
hold on
plot(s,diff(2,:),"b",'linewidth',2)
plot(s,diff(3,:),"g",'linewidth',2)
xlabel("s")
ylabel("\Delta")
%legend("\Delta x", "\Delta y", "\Delta z")
hold off
%}

%{
% FIGURE 6 IN PAPER
figure(4)
subplot(2,1,1)
plot(s,magn(1,:),'r','linewidth',2)
hold on
plot(s,magn(3,:),'b','linewidth',2)
plot(s,magn(4,:),'g','linewidth',2)
xlabel("s")
ylabel("|X|")
%ylim([0.75 5.75])
title("Expanded Method")
%legend("Tangent Vector","Principal Normal Vector", "Binormal Vector","Location","Northwest")
hold off
subplot(2,1,2)
plot(s,magn(2,:),'r','linewidth',2)
hold on
plot(s,magn(3,:),'b','linewidth',2)
plot(s,magn(5,:),'g','linewidth',2)
xlabel("s")
ylabel("|X|")
%ylim([0.75 5.75])
title("Non-Expanded Method")
%legend("Tangent Vector","Principal Normal Vector", "Binormal Vector","Location","Northwest")
hold off
%

%
% FIGURE 7 IN PAPER
figure(5)
subplot(2,1,1)
plot(s,dot_tnb(1,:),'r','linewidth',2)
hold on
plot(s,dot_tnb(2,:),'g','linewidth',2)
plot(s,dot_tnb(3,:),'b','linewidth',2)
xlabel("s")
ylabel("X_{1}\cdot X_{2}")
%ylim([-1 1.2])
title("Expanded Method")
%legend("Tangent Vector","Principal Normal Vector", "Binormal Vector","Location","Northwest")
hold off
subplot(2,1,2)
plot(s,dot_tnb(4,:),'r','linewidth',2)
hold on
plot(s,dot_tnb(5,:),'g','linewidth',2)
plot(s,dot_tnb(6,:),'b','linewidth',2)
xlabel("s")
ylabel("X_{1}\cdot X_{2}")
%ylim([-1 1.2])
title("Non-Expanded Method")
%legend("Tangent Vector","Principal Normal Vector", "Binormal Vector","Location","Northwest")
hold off
%}


%{
% PLOTTING F_cm FOR EXPANDED AND NON-EXPANDED METHODS
figure(6)
subplot(1,3,1)
plot(s,fun_curve_0_exp, "r",'linewidth',2) % f_c0 function obtained from expanded method
hold on
plot(s,fun_curve_0_nonexp, "b",'linewidth',2) % f_c0 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C0}(s)")
%ylim([-.05 7])
%xlim([1 5])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(1,3,2)
plot(s,fun_curve_1_exp, "r",'linewidth',2) % f_c1 function obtained from expanded method
hold on
plot(s,fun_curve_1_nonexp, "b",'linewidth',2) % f_c1 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C1}(s)")
%ylim([-.05 4.5])
%xlim([1 5])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(1,3,3)
plot(s,fun_curve_2_exp, "r",'linewidth',2) % f_c2 function obtained from expanded method
hold on
plot(s,fun_curve_2_nonexp, "b",'linewidth',2) % f_c2 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C2}(s)")
%ylim([-.05 7])
%xlim([1 5])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
%}


%{
% PLOTTING F_tm FOR EXPANDED AND NON-EXPANDED METHODS
figure(7)
subplot(1,3,1)
plot(s,fun_tang_exp_0, "r",'linewidth',2) % f_t0 function obtained from expanded method
hold on
plot(s,fun_tang_nonexp_0, "b",'linewidth',2) % f_t0 function obtained from non-expanded method
xlabel("s")
ylabel("f_{T0}(s)")
%xlim([1 5])
%ylim([-.05 2])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(1,3,2)
plot(s,fun_tang_exp_1, "r",'linewidth',2) % f_t1 function obtained from expanded method
hold on
plot(s,fun_tang_nonexp_1, "b",'linewidth',2) % f_t1 function obtained from non-expanded method
xlabel("s")
ylabel("f_{T1}(s)")
%xlim([1 5])
%ylim([-.05 1])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(1,3,3)
plot(s,fun_tang_exp_2, "r",'linewidth',2) % f_t2 function obtained from expanded method
hold on
plot(s,fun_tang_nonexp_2, "b",'linewidth',2) % f_t2 function obtained from non-expanded method
xlabel("s")
ylabel("f_{T2}(s)")
%xlim([1 5])
%ylim([-.05 1])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
%}



%{
% PLOTTING F_cm FOR EXPANDED AND NON-EXPANDED METHODS
figure(8)
subplot(2,3,1)
plot(s,fun_tang_exp_0, "r",'linewidth',2) % f_t0 function obtained from expanded method
hold on
plot(s,fun_tang_nonexp_0, "b",'linewidth',2) % f_t0 function obtained from non-expanded method
xlabel("s")
ylabel("f_{T0}(s)")
%xlim([1 5])
%ylim([-.05 2])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(2,3,2)
plot(s,fun_tang_exp_1, "r",'linewidth',2) % f_t1 function obtained from expanded method
hold on
plot(s,fun_tang_nonexp_1, "b",'linewidth',2) % f_t1 function obtained from non-expanded method
xlabel("s")
ylabel("f_{T1}(s)")
%xlim([1 5])
%ylim([-.05 1])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(2,3,3)
plot(s,fun_tang_exp_2, "r",'linewidth',2) % f_t2 function obtained from expanded method
hold on
plot(s,fun_tang_nonexp_2, "b",'linewidth',2) % f_t2 function obtained from non-expanded method
xlabel("s")
ylabel("f_{T2}(s)")
%xlim([1 5])
%ylim([-.05 1])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(2,3,4)
plot(s,fun_curve_0_exp, "r",'linewidth',2) % f_c0 function obtained from expanded method
hold on
plot(s,fun_curve_0_nonexp, "b",'linewidth',2) % f_c0 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C0}(s)")
%ylim([-.05 7])
%xlim([1 5])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(2,3,5)
plot(s,fun_curve_1_exp, "r",'linewidth',2) % f_c1 function obtained from expanded method
hold on
plot(s,fun_curve_1_nonexp, "b",'linewidth',2) % f_c1 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C1}(s)")
%ylim([-.05 4.5])
%xlim([1 5])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off
subplot(2,3,6)
plot(s,fun_curve_2_exp, "r",'linewidth',2) % f_c2 function obtained from expanded method
hold on
plot(s,fun_curve_2_nonexp, "b",'linewidth',2) % f_c2 function obtained from non-expanded method
xlabel("s")
ylabel("f_{C2}(s)")
%ylim([-.05 7])
%xlim([1 5])
%legend("Expanded method","Non-Expanded method",'Location','northwest')
hold off


%}
