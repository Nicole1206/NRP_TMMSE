function [Rglc_mwf,alpha_mwf, belta_mwf] = MWF_TMMSE_function(x,R,Uxx, Sxx,Ext)

M = size(x, 1);
N = size(x, 2);

% Rnn = mean(diag(Sxx(Ext + 1:end, Ext + 1:end))) * eye(M);
Rnn = mean(diag(Sxx(end-Ext:end, end-Ext:end))) * eye(M);

% Sxx_sig = zeros(M, M);
% Snn = zeros(M, M);
% 
% Sxx_sig(1:Ext, 1:Ext) = Sxx(1:Ext, 1:Ext);
% Snn(end-Ext:end, end-Ext:end) = Sxx(end-Ext:end, end-Ext:end);
% Rxx = Uxx * Sxx_sig * Uxx';
% Rnn = Uxx * Snn * Uxx';

Rxx_temp = R - Rnn;
[V, D] = eig(Rxx_temp);
realDiag = real(diag(D));
imagDiag = imag(diag(D));
realDiag = max(realDiag, 0.0);
D = diag(realDiag + 1j*imagDiag);
Rxx = V * D * V';
mwf = inv(R) * (Rxx); % standrad MWF
% mwf = inv(R) * Rxx;
% u = 1;
% mwf = (Rxx + u*Rnn).^(-1) * Rxx; % SDW-MWF
% mwf = inv(Rnn) * Rxx * 1 / (u + trace( inv(Rnn) * Rxx )); % Rank-one SDW-MWF
% r0 = Rxx(1,1);
% mwf = inv(Rnn) * Rxx * r0 / (u*r0 + trace( inv(Rnn) * (Rxx(:,1)*Rxx(1,:)) )); % Saptial Prediction SDW-MWF

y = mwf' * x;


Ryy = zeros(M,M);
for l = 1 : size(y,2)
    Ryy = Ryy + y(:,l) * y(:,1)';
end
Ryy = Ryy / size(y,2);


[Uyy,Syy] = svd(Ryy);
Rpri = Uyy * diag([ones(1,Ext) zeros(1,M-Ext)]) * Uyy';
[Rglc_mwf,alpha_mwf, belta_mwf] = GLCfunction_rec(Ryy,Rpri,R,M,N,Ext,y);