% ------------------------------------------------------------------- 
% SVD-based mixed-type Kalman Filter 
%         Method: SVD-based implementation called V-Lambda filter
%           Type: Covariance filtering
%      Recursion: Riccati recursion
%           Form: Two stages, a posteriori form
%        Authors: Oshman Y., Bar-Itzhack, I.Y. (1986)
% Implementation: Maria Kulikova 
% ------------------------------------------------------------------- 
% References: (see Table 2, a posteriori gain calculation)
%   1. Oshman Y., Bar-Itzhack, I.Y. (1986)
%      Square root filtering via covariance and information eigenfactors.
%      Automatica, 22(5):599-604. DOI: https://doi.org/10.1016/0005-1098(86)90070-1 
% ------------------------------------------------------------------- 
% Input:
%     matrices        - system matrices F,H,Q etc
%     initials_filter - initials x0,P0
%     measurements    - measurements (where y(t_k) is the k-th column)
% Output:
%     neg_LLF     - negative log LF
%     hatX        - filtered estimate (history) 
%     hatDP       - diag of the filtered error covariance (history)
% ------------------------------------------------------------------- 

function [neg_LLF,hatX,hatDP] = Riccati_KF_VLambda1(matrices,initials_filter,measurements)
   [F,G,Qsys,H,R] = deal(matrices{:});         % get system matrices
            [X,P] = deal(initials_filter{:});  % get initials for the filter 
   
        [m,n]  = size(H);                % get dimensions
       N_total = size(measurements,2);   % number of measurements
          hatX = zeros(n,N_total+1);     % prelocate for efficiency
         hatDP = zeros(n,N_total+1);     % prelocate for efficiency

     [~,DQ,QQ] = svd(Qsys);      % SVD for process noise covariance
        GQsqrt = G*QQ*DQ.^(1/2); % to compute once the new G*Q^{1/2}
     [~,DR,RQ] = svd(R);         % SVD for measurement noise covariance
        R_sqrt = RQ*DR.^(1/2); clear DR; % to compute once
     [~,DP,QP] = svd(P); DPsqrt = DP.^(1/2); clear DP;

   neg_LLF = 1/2*m*log(2*pi)*N_total;                 % initial value for the neg Log LF
hatX(:,1)  = X; hatDP(:,1) = diag(QP*(DPsqrt^2)*QP'); % save initials at the first entry
for k = 1:N_total                       
   [X,DPsqrt,QP]        = svd_predict(X,DPsqrt,QP,F,GQsqrt); 
   [X,DPsqrt,QP,ek,Rek] = svd_update(X,DPsqrt,QP,measurements(:,k),H,R_sqrt);
   
   neg_LLF = neg_LLF+1/2*log((det(Rek)))+1/2*ek'/Rek*ek;
   hatX(:,k+1)= X; hatDP(:,k+1) = diag(QP*DPsqrt^(2)*QP'); % save estimates  
 end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Time update: a priori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,DPsqrt,QP] = svd_predict(X,DPsqrt,QP,F,GQsqrt)
      PreArray   = [F*QP*DPsqrt, GQsqrt];
           [m,n] = size(PreArray);        
     if n<m, error('The matrix is fat, i.e. m>n'); end;

     [QP,DPsqrt] = svd(PreArray,'econ'); 
               X = F*X;                  
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Measurement update: a posteriori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,DPsqrt,QP,ek,Rek] = svd_update(X,DPsqrt,QP,z,H,R_sqrt)
     Rek = R_sqrt*R_sqrt' + H*(QP*(DPsqrt^2)*QP')*H';  % residual covariance for LLF
     %_______________________________________________
     PreArray   = [QP/DPsqrt, H'/R_sqrt'];
          [m,n] = size(PreArray);        
     if n<m, error('The matrix is fat, i.e. m>n'); end;
     [QP,invDPsqrt] = svd(PreArray,'econ'); 
         DPsqrt = inv(invDPsqrt);                   
 
     Gain = (QP*DPsqrt*DPsqrt*QP')*H'/R_sqrt'/R_sqrt;   % Kalman Gain
       ek = z - H*X;                                    % residual for LLF
        X = X + Gain*ek;                                % Filtered state estimate        
end
