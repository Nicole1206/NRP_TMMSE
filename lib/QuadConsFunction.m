function [Wqc, DLqc] = QuadConsFunction(Ang_left, Ang_right, TestPoint, kd, M, R, DLqc)


AngTest = Ang_left : (Ang_right-Ang_left)/TestPoint : Ang_right;
adS_tp = [(exp(1j * kd * sin(Ang_left) * (1 : M)')), (exp(1j * kd * sin(Ang_right) * (1 : M)'))]; % The two point quadratic SV

%loop
% DLqc = 1; % initiate the DLL
StepSize = 2;

while 1
    Rqc = R + DLqc*eye(M);
    %             Rqc = eye(M);
    R_temp = inv( adS_tp' * inv(Rqc) * adS_tp );
    r0 = R_temp(1, 1);
    r1 = R_temp(2, 2);
    r2 = abs(R_temp(1, 2));
    
    if real(r2) >= 0
        Phrase_beta = angle(R_temp(1, 2));
    else
        Phrase_beta = angle(R_temp(1, 2)) + pi;
    end
    Phrase_fai = pi - Phrase_beta;
    
    % KKT
    if r2/r0 <= 1
        rho0 = 1;
    else
        rho0 = r2/r0;
    end
    
    if r2/r1 <= 1
        rho1 = 1;
    else
        rho1 = r2/r1;
    end
    
    Wqc = inv(Rqc) * adS_tp * R_temp * [rho0, rho1*exp(1j*Phrase_fai)].';
    
    EndFlag = 0;
    for AngIdx = 2 : length(AngTest) - 1
        adS_med = exp(1j * kd * sin(AngTest(AngIdx)) * (1 : M)');
        if abs(adS_med'*Wqc) >= 1
            EndFlag = EndFlag + 1;
            %                     AngTest = AngTest(2:end);
        else
            DLqc = StepSize * DLqc;
            break;
        end
    end
    
    if EndFlag == length(AngTest) - 2
        break;
    end
    
end