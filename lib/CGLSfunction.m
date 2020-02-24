function  [Wcg, iter, rs] = CGLSfunction(ConsMat, Gv, x, tol_dB, DLcg, SCtype) 
% CGLSfunction: Conjugate Gradient Least Square Method.
 %   ConsMat is the constraint matrix of the LCMV beamformer; 
 %   Gv is the Gain vector.
 %   tol_dB is the threshold for stopping criteria.
 %   DLcg is the diagonal loading value for beamforming.
 %   Stopping criteria type: GDP for generalized discrepancy principle; 
 %                           EE for error estimation; 
 %                           RVE for Ritz value estimation.
 
 
N = size(x,2);
M = size(x, 1);
BlockMat = [eye(M-size(ConsMat, 2)); zeros(1, M-size(ConsMat, 2))] - ConsMat*(ConsMat'*ConsMat)^(-1)*ConsMat'*[eye(M-size(ConsMat, 2)); zeros(1, M-size(ConsMat, 2))];
W_qui = ConsMat*(ConsMat'*ConsMat)^(-1)*Gv;

d_sig = W_qui' * x;
X_bs = BlockMat' * x;

Rx_bs = 1/N * (X_bs) * X_bs';   % the covariance matrix of the blocked signal
r_x1x2 = 1/N * X_bs * d_sig'; % the cross-correlation vector between x1 and d_sig

% CGLS iteration
% initiation
tol = 10^(tol_dB/20);
W_gsc = zeros(M-size(ConsMat, 2), 1); % the adaptive weight vector of the blocking matrix
% res_LS = d_sig' - X_bs' * W_gsc;
% norm_LSres = norm(res_LS);
% r(1) = norm_LSres;

res_NE = r_x1x2 - Rx_bs * W_gsc;
norm_NEres = norm(res_NE);
rs(1) = norm(res_NE);

dir_search = res_NE;
squa_norm_NEres = (norm(res_NE))^2;
iter = 0;
Q = 0;
P = 0;

% loop
if strcmp(SCtype, 'GDP')
    while norm_NEres > tol
        iter = iter + 1;
%         vec_temp = X_bs' * dir_search;
        step_alpha = squa_norm_NEres / (dir_search' * (Rx_bs + DLcg*eye(M-size(ConsMat, 2))) * dir_search);

        W_gsc = W_gsc + step_alpha * dir_search;
%         res_LS = res_LS - step_alpha * vec_temp;
%         norm_LSres = norm(res_LS);
%         r(iter + 1) = norm_LSres;
        res_NE = res_NE - step_alpha * (Rx_bs + DLcg*eye(M-size(ConsMat, 2))) * dir_search;
        norm_NEres = norm(res_NE);
        rs(iter + 1) = norm_NEres;
        
        squa_norm_NEres_last = squa_norm_NEres;
        squa_norm_NEres = (norm(res_NE))^2;
        step_beta = squa_norm_NEres / squa_norm_NEres_last;
        dir_search = res_NE + step_beta * dir_search;
    end
    
    
elseif strcmp(SCtype, 'EE')
    disp('warning: I have not finished the code of EE')
    while 1
        iter = iter + 1;
        vec_temp = X_bs' * dir_search;
        step_alpha = squa_norm_NEres / (vec_temp'*vec_temp  + DLcg * (dir_search)'*dir_search);
        W_gsc = W_gsc + step_alpha * dir_search;
        
        squa_norm_LSres = (norm(res_LS))^2;
        res_LS = res_LS - step_alpha * vec_temp;
        res_temp = res_temp + DLcg * step_alpha * dir_search;
        res_NE = X_bs * res_LS - res_temp;
        
        squa_norm_NEres_last = squa_norm_NEres;
        squa_norm_NEres = (norm(res_NE))^2;
        step_beta = squa_norm_NEres / squa_norm_NEres_last;
        dir_search = res_NE + step_beta * dir_search;
        
        % Stopping Criteria
        Q = 1 + step_beta * Q;
        P = P + step_alpha * Q;
        squa_norm_LSres = squa_norm_LSres - step_alpha * squa_norm_NEres_last;
        if sqrt(P * squa_norm_LSres) <=  tol
            break;
        end
    end
    
    
elseif strcmp(SCtype, 'RVE')
    disp('warning: I have not finished the code of RVE')
    while 1
        iter = iter + 1;
        vec_temp = X_bs' * dir_search;
        step_alpha = squa_norm_NEres / (vec_temp'*vec_temp  + DLcg * (dir_search)'*dir_search);
        W_gsc = W_gsc + step_alpha * dir_search;
        res_LS = res_LS - step_alpha * vec_temp;
        res_temp = res_temp + DLcg * step_alpha * dir_search;
        res_NE = X_bs * res_LS - res_temp;
        
        squa_norm_NEres_last = squa_norm_NEres;
        squa_norm_NEres = (norm(res_NE))^2;
        step_beta = squa_norm_NEres / squa_norm_NEres_last;
        dir_search = res_NE + step_beta * dir_search;
        
        % Stopping Criteria
        Q = 1 + step_beta * Q;
        P = P + step_alpha * Q;
        if step_alpha * Q <= 1/ tol
            break;
        end
    end
else
    disp('Error: Input the wrong stopping criteria')
end


Wcg = W_qui - BlockMat * W_gsc;

end