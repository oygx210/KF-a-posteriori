% ------------------------------------------------------------------- 
% Classical Kalman Filter 
%         Method: Conventional implementation
%           Type: Covariance filtering
%      Recursion: Riccati recursion
%           Form: Two stages, a posteriori form
%        Authors: R.E. Kalman 
% Implementation: Maria Kulikova 
% ------------------------------------------------------------------- 
% References:
%   Kalman R.E. 
%   A new approach to linear filtering and prediction problems. 
%   Journal of basic Engineering. 1960 Mar, 82(1):35-45.
%
%   See also implementation in
%   Grewal, M.S., Andrews, A.P. (2015) 
%   Kalman filtering: theory and practice using MATLAB 
%   Prentice-Hall, New Jersey, 4th edn. 
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

function [neg_LLF,hatX,hatDP] = Riccati_KF_standard(matrices,initials_filter,measurements)
   [F,G,Q,H,R] = deal(matrices{:});         % get system matrices
         [X,P] = deal(initials_filter{:});  % get initials for the filter 
   
        [m,n]  = size(H);                   % dimensions
       N_total = size(measurements,2);      % number of measurements
         hatX  = zeros(n,N_total+1);        % prelocate for efficiency
         hatDP = zeros(n,N_total+1);        % prelocate for efficiency
       neg_LLF = 1/2*m*log(2*pi)*N_total;   % set initial value for the neg Log LF
 hatX(:,1)  = X; hatDP(:,1) = diag(P);      % save initials at the first entry
 for k = 1:N_total                
      [X,P]        = kf_predict(X,P,F,G,Q); 
      [X,P,ek,Rek] = kf_update(X,P,measurements(:,k),H,R);     
      
      neg_LLF = neg_LLF+1/2*log((det(Rek)))+1/2*ek'/Rek*ek;
      hatX(:,k+1)  = X; hatDP(:,k+1) = diag(P); 
  end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Time update: a priori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,P] = kf_predict(X,P,F,G,Q)
    X = F*X;                  % Predicted state estimate  
    P = F*P*F' + G*Q*G';      % Predicted error covariance 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Measurement update: a posteriori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,P,residual,cov_residual] = kf_update(X,P,z,H,R)

  residual     = z - H*X;                 % residual
  cov_residual = R + H*P*H';              % residual covariance 
  Kalman_gain  = P*H'/cov_residual;       % Filter gain

  X = X + Kalman_gain*residual;           % Filtered state estimate
  P = P - Kalman_gain*H*P;                % Filtered error covariance
end
