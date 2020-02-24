function Inter = InterGenerateFunction(M,N,Ener,adS)

Inter = zeros(M,N);
% x = zeros(1,length(N));
% for l = 1 : N,
%     x(l) = wgn(1,1,Ener,'complex');
% end;
x = wgn(1,N,Ener,'complex');
% x = wgn(1,N,Ener,'real');


Inter = adS * x;