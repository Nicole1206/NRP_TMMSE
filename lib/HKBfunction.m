function [rhoLS, Rhkb] = HKBfunction(ad, R, M)

        [B,r] = qr(ad);
%         [B,v] = eig(eye(M) - ad * ad'/M);
        B = B(: , 2 : end);
        X = sqrtm(R) * B;
        b = sqrtm(R) * ad / M;
        eta = (X'*X)^(-1)*X'*b;
        deltaLS = norm(X * eta - b);
        rhoLS = (M - 1)*deltaLS^2 / (norm(eta))^2;
        Rhkb = R + rhoLS * eye(M);
end