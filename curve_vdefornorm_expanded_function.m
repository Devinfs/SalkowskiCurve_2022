function [curve, tang, norm, bino]=curve_vdefornorm_expanded_function(s, kap)
% This functions goes with the fourth draft of the paper.


% initial s
s0=1;
% t interval
t=kap*log(s);

% B matrix (array of coefficients of F_Nm)(EQ. 4.2)
B=[[1,0,0, -1/(factorial(3)*kap),-1/(factorial(4)*(kap^2)),(1/factorial(5))*((-13/(kap^3))+(1/kap))];[0,1,0, -2/factorial(3),-5/(factorial(4)*kap),(1/factorial(5))*((-17/(kap^2))+4)];[0,0,1,2/(factorial(3)*kap),(1/factorial(4))*((2/(kap^2))-4),(1/factorial(5))*((2/(kap^3))-(18/kap))]];

% functions needed for normal vector (EQ. 4.3)
F_N0= B(1,1) + B(1,2)*t + B(1,3)*t.^2 + B(1,4)*t.^3 + B(1,5)*t.^4 + B(1,6)*t.^5;
F_N1= B(2,1) + B(2,2)*t + B(2,3)*t.^2 + B(2,4)*t.^3 + B(2,5)*t.^4 + B(2,6)*t.^5;
F_N2= B(3,1) + B(3,2)*t + B(3,3)*t.^2 + B(3,4)*t.^3 + B(3,5)*t.^4 + B(3,6)*t.^5;

% PRINCIPAL NORMAL VECTOR (THM. 3)
norm(1,:)=-F_N1 - ((1/(2*kap))*F_N2);
norm(2,:)=F_N0 - F_N2;
norm(3,:)=F_N1;

% G_mp matrix (array of coeffecients of F_Tm) (EQ. 5.2)
G_mp= [[1, 1/(2*kap), 1/(factorial(3)*kap^2), (1/(factorial(4)))*((1/kap^3)-(1/kap)), (1/factorial(5))*((1/kap^4)-(9/kap^2)),0,0,0,0,0]; [0, 1/2, 1/(3*kap), ((1/(8*kap^2))-(1/12)),((1/(30*kap^3))-(13/(120*kap^2))),0,0,0,0,0];[0,0,1/3,1/(3*kap),((-1/30)+(11/(60*kap^2))), ((1/(45*kap^3))-(19/(360*kap)))],0,0,0,0];

% functions needed for tangent vector (EQ. 5.5)
F_T0= G_mp(1,1)*t + G_mp(1,2)*t.^2 + + G_mp(1,3)*t.^3 + G_mp(1,4)*t.^4 + G_mp(1,5)*t.^5 + G_mp(1,6)*t.^6;
F_T1= G_mp(2,1)*t + G_mp(2,2)*t.^2 + + G_mp(2,3)*t.^3 + G_mp(2,4)*t.^4 + G_mp(2,5)*t.^5 + G_mp(2,6)*t.^6;
F_T2= G_mp(3,1)*t + G_mp(3,2)*t.^2 + + G_mp(3,3)*t.^3 + G_mp(3,4)*t.^4 + G_mp(3,5)*t.^5 + G_mp(3,6)*t.^6;

% TANGENT VECTOR (THM. 3)
tang(1,:)=1-F_T1 - (1/(2*kap))*F_T2;
tang(2,:)=F_T0 - F_T2;
tang(3,:)=F_T1;

% BINORMAL VECTOR (THM. 3)
bino(1,:)=tang(2,:).*norm(3,:)-tang(3,:).*norm(2,:);
bino(2,:)=tang(3,:).*norm(1,:)-tang(1,:).*norm(3,:);
bino(3,:)=tang(1,:).*norm(2,:)-tang(2,:).*norm(1,:);


% ET matrix (array with values needed to calculate G_mp values) (PG. 18)
ET_exp=[[1, 1, 1, 1, 1, 1]; [1, 1, 1, 1, 1, 1]; [0, 1, 1, 1, 1, 1]; [0, 0, 1, 1, 1, 1]; [0, 0, 0, 1, 1, 1]; [0, 0, 0, 0, 1, 1]; [0, 0, 0, 0, 0, 1]];

% H_mq matrix (APPENDIX B, PG. 18)
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

% The next two "%" relate to the final functions needed for curve 

% integration constants of main functions (EQ. 6.6 FOR THM. 3)
F_curve0_cons= s0.*(H_mq(1,1)+ H_mq(1,2)*log(s0)+ H_mq(1,3)*log(s0).^2 + H_mq(1,4)*log(s0).^3 + H_mq(1,5)*log(s0).^4 + + H_mq(1,6)*log(s0).^5 + H_mq(1,7)*log(s0).^6);
F_curve1_cons= s0.*(H_mq(2,1)+ H_mq(2,2)*log(s0)+ H_mq(2,3)*log(s0).^2 + H_mq(2,4)*log(s0).^3 + H_mq(2,5)*log(s0).^4 + + H_mq(2,6)*log(s0).^5 + H_mq(2,7)*log(s0).^6);
F_curve2_cons= s0.*(H_mq(3,1)+ H_mq(3,2)*log(s0)+ H_mq(3,3)*(log(s0).^2) + H_mq(3,4)*log(s0).^3 + H_mq(3,5)*log(s0).^4 + + H_mq(3,6)*log(s0).^5 + H_mq(3,7)*log(s0).^6);

% main functions (EQ. 6.6 FOR THM. 3)
F_curve0= s.*(H_mq(1,1)+ H_mq(1,2)*log(s)+ H_mq(1,3)*log(s).^2 + H_mq(1,4)*log(s).^3 + H_mq(1,5)*log(s).^4 + + H_mq(1,6)*log(s).^5 + H_mq(1,7)*log(s).^6) - F_curve0_cons;
F_curve1= s.*(H_mq(2,1)+ H_mq(2,2)*log(s)+ H_mq(2,3)*log(s).^2 + H_mq(2,4)*log(s).^3 + H_mq(2,5)*log(s).^4 + + H_mq(2,6)*log(s).^5 + H_mq(2,7)*log(s).^6) -F_curve1_cons;
F_curve2= s.*(H_mq(3,1)+ H_mq(3,2)*log(s)+ H_mq(3,3)*(log(s).^2) + H_mq(3,4)*log(s).^3 + H_mq(3,5)*log(s).^4 + + H_mq(3,6)*log(s).^5 + H_mq(3,7)*log(s).^6) -F_curve2_cons;

% curve (THM. 3)
curve(1,:)= s-s0-(1/(2*kap))*F_curve2 - F_curve1; %COEFFICIENT OF VECTOR A
curve(2,:)=F_curve0-F_curve2; %COEFFICIENT OF VECTOR B
curve(3,:)=-F_curve1; %COEFFICIENT OF VECTOR C

end