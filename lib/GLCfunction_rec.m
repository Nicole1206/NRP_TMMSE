function [Rglc_new,alpha0,belta0, u, rho, gamma] = GLCfunction_rec(Ryy,Rpri,R,M,N,L,y)


u = trace(Rpri'* Ryy) / L;
quadruple = 0;
for l = 1 : N,
    quadruple = quadruple + norm(y(:,l))^4;
end;
rho = 1/N^2 * quadruple - 1/N * trace(Ryy' * Ryy);
diff = u * Rpri - Ryy;
gamma = trace(diff' * diff);
belta0 = 1 - rho / gamma;
alpha0 = u*(1-belta0);

% v = trace(R) / M;
% gamma = sum(sum(abs(v * eye(M) - R).^2));
% rho = 1/N^2 * sum(sum(abs(x),2).^4) - 1/N * sum(sum(abs(R).^2));
% rho1 = sum(sum(abs(x - mean(x,2) * ones(1,N)).^2));
% rho = 0;
% for l = 1 : N,
%     rho = rho + 1/N^2 * sum(sum(abs(x(:,l) * x(:,l)' - R)).^2);
% end;
% rho = 0;
% for l = 1 : N,
%     rho = rho + 1/N^2 * (x(:,l)' * x(:,l))^2;
% end;
% rho = rho - 1/N * sum(sum(abs(R).^2));
% 
% alpha0 = min(v*rho/gamma,v);
% belta0 = 1 - alpha0 / v;
% 
% R = alpha0 * eye(M) + belta0 * R;
Rglc_new = alpha0 / belta0 * eye(M) + R;