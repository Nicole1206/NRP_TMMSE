function [Rglc_LS,alpha_LS, belta_LS, u, rho, gamma, Rpri_LS] = LS_TMMSE_function(Hm,Hl,x,R,Ext)

M = size(x, 1);
N = size(x, 2);

yls = zeros(M, N);
for ll = 1 : M
    x_single = x(ll,:);
    x_hankel = hankel(x_single(1 : Hl), x_single(Hl : end));
    [U,S,V] = svd(x_hankel);
    U1 = U(:,1:Ext);
    V1 = V(:,1:Ext);
    S1 = S(1:Ext,1:Ext);
    y_hankel = U1*S1*V1';
    y_single = zeros(1, N);
    for Idx = 1 : N
        A = max(1,Idx - Hl + 1);
        B = min(Hm, Idx);
        sum_hankel = 0;
        for kIdx = A : B
            sum_hankel = sum_hankel + y_hankel(Idx - kIdx + 1, kIdx );
        end
        y_single(Idx) = sum_hankel/(B-A+1);
    end
    yls(ll,:) = y_single;
end

%             yls = yls';
Ryy_LS = zeros(M,M);
for l = 1 : size(yls,2)
    Ryy_LS = Ryy_LS + yls(:,l) * yls(:,l)';
end
Ryy_LS = Ryy_LS / size(yls,2);


[Uyy_LS,Syy_LS] = svd(Ryy_LS);
Rpri_LS = Uyy_LS * diag([ones(1,Ext) zeros(1,M-Ext)]) * Uyy_LS';
[Rglc_LS,alpha_LS, belta_LS, u, rho, gamma] = GLCfunction_rec(Ryy_LS,Rpri_LS,R,M,N,Ext,yls);