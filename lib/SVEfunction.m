function ad_est = SVEfunction(Ang_left,Ang_right,Ang_interval,kd, R)

M = size(R,1);

% the complement of the target sector
Ang_comp = [-pi : Ang_interval : Ang_left, Ang_right : Ang_interval : pi];
R_Ang_sv = zeros(M, M);
for AngIdx = 1 : length(Ang_comp)
    Ang_sv = exp(1j * kd * sin(Ang_comp(AngIdx)) * (1 : M)');
    R_Ang_sv = R_Ang_sv + Ang_sv * Ang_sv';
end
Sector_comp = R_Ang_sv * Ang_interval;

Ang_tar = Ang_left : Ang_interval : Ang_right;
th_array = zeros(1, length(Ang_tar));
for AngIdx = 1 : length(Ang_tar)
    Ang_sv = exp(1j * kd * sin(Ang_tar(AngIdx)) * (1 : M)');
    th_array(AngIdx) = Ang_sv' * Sector_comp * Ang_sv;
end
th = max(th_array);

%         % value of the term ad' * Sector_comp * ad for different angles
%         Ang_all = -pi : Ang_interval : pi;
%         th_all_array = zeros(1, length(Ang_all));
%         for AngIdx = 1 : length(Ang_all)
%             Ang_sv = exp(1j * kd * sin(Ang_all(AngIdx)) * (1 : M)');
%             th_all_array(AngIdx) = Ang_sv' * Sector_comp * Ang_sv;
%         end
%         plot(Ang_all/pi*180, abs(th_all_array));

% solve the SDP
cvx_begin
variable A(M,M) complex semidefinite;
minimize( abs(trace(R\A)) );
subject to
trace(A) == M;
abs(trace(Sector_comp*A)) <= abs(th);
cvx_end

% estimate the steering vector
[Usve, Ssve] = svd(A);
ad_est = Usve(:,1);
