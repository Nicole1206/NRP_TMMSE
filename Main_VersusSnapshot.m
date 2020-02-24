clear all;
close all;
clc;

addpath(genpath('lib'));

c = 344;
f = 1000;
d = c/f/2; %half-wavelength
kd = 2 * pi * f/c * d;

M = 20;

ulaPos = zeros(M, 3);
ulaPos(:, 2) = (1 : M) * d;

N_list = [10 : 5 : 20, 30: 20 : 70, 100 : 50 :200, 300 : 100 : 500, 1000 : 1000 : 4000];
% N = 1000;
flag = 1; % Perturbation existing
% delta = 0 / 180 * pi;
% delta_list = [0:0.1:3]/180*pi;
% SNR_list = -10:5:40;


for NIdx = 1 : length(N_list)
    N = N_list(NIdx);
    NIdx
    for m = 1 : 1,
        
        EnerW = 0;
        EnerS = 10;
        EnerI = 20;
        % EnerI2 = 20;
        
        Int = 3;
        DOAS = 0/180*pi;
        DOAI = (30 + (0:(Int - 1))*15)/180*pi;
        
        % position perturbation in X-Y plane
        if flag == 1
            posPtb = 0.05 * rand * (2*d);
            theta = pi * (2*rand(M, 1) - 1);
            ulaPosPtb = ulaPos;
            ulaPosPtb(:,1) = ulaPosPtb(:,1) + posPtb * cos(theta);
            ulaPosPtb(:,2) = ulaPosPtb(:,2) + posPtb * sin(theta);
        else
            ulaPosPtb = ulaPos;
        end
        
        delta = 5 * (2*rand - 1) * pi / 180;
%         adS = exp(1j * kd * sin(DOAS) * (1 : M)');
        adS = SteerVectorGenerate(ulaPosPtb, DOAS + delta);
        Signal = InterGenerateFunction(M,N,EnerS,adS);
        
        adI = zeros(M,1,Int);
        Inter = zeros(M,N,Int);
        for IIdx = 1 : Int
            adI(:,:,IIdx) = SteerVectorGenerate(ulaPosPtb, DOAI(IIdx));
            Inter(:,:,IIdx) = InterGenerateFunction(M,N,EnerI,adI(:,:,IIdx));
        end
        
        WGN = wgn(M,N,EnerW,'complex');
        
        x = Signal + sum(Inter,3) + WGN;
        Q = sum(10^(EnerI/10) * adI(:,:,IIdx) * (adI(:,:,IIdx))',3)  + 10^(EnerW/10) * eye(M);
        R0 = 10^(EnerS/10) * (adS) * (adS)' + Q;
        
        ad = SteerVectorGenerate(ulaPos, DOAS);
        [SINR_opt(m),SOI_opt(m)] = PerformanceF(EnerS,Q,R0,R0,1,adS,adS);
        
        R = zeros(M,M);
        for l = 1 : size(x,2),
            R = R + x(:,l) * x(:,l)';
        end;
        R = R / size(x,2);
        [Uxx,Sxx] = svd(R);
        [SINR_SCB(m),SOI_SCB(m)] = PerformanceF(EnerS,Q,R,R,1,ad,adS);
        
        %method of SVE
        Ang_left = DOAS - 5*pi/180;
        Ang_right = DOAS + 5*pi/180;
        Ang_interval = 0.1 * pi / 180;       
        ad_est = SVEfunction(Ang_left,Ang_right,Ang_interval,kd, R);
        
       %% LS method
        Ext = Int + 1;
        Hm = 2 * 3;
        Hl = N + 1 - Hm;        
        [Rglc_LS,alpha_LS, belta_LS] = LS_TMMSE_function(Hm,Hl,x,R,Ext);
        
        DL_glc_ls(m) = abs(alpha_LS/belta_LS);  
        [SINR_GLC_LS(m),SOI_GLC_LS(m)] = PerformanceF(EnerS,Q,Rglc_LS,R,1,ad,adS);
        
        % the SVE of LS-TMMSE      
%         ad_ls_sve = SVEfunction(Ang_left,Ang_right,Ang_interval,kd, Rglc_LS);
        [SINR_GLC_LS_SVE(m),SOI_GLC_LS_SVE(m)] = PerformanceF(EnerS,Q,Rglc_LS,R,1,ad_est,adS);
        
        %% method of Winer Filter   
        [Rglc_mwf,alpha_mwf, belta_mwf] = MWF_TMMSE_function(x,R,Uxx,Sxx,Ext);
       
        DL_glc_mwf(m) = abs(alpha_mwf/belta_mwf);
        [SINR_GLC_MWF(m),SOI_GLC_MWF(m)] = PerformanceF(EnerS,Q,Rglc_mwf,R,1,ad,adS);
        
%         ad_mwf_sve = SVEfunction(Ang_left,Ang_right,Ang_interval,kd, Rglc_mwf);
        [SINR_GLC_MWF_SVE(m),SOI_GLC_MWF_SVE(m)] = PerformanceF(EnerS,Q,Rglc_mwf,R,1,ad_est,adS);
        
        %% method of GLC
        Ext = M;
        Rpri = Uxx * diag([ones(1,Ext) zeros(1,M-Ext)]) * Uxx';
        [Rglc,alpha0,belta0] = GLCfunction_rec(R,Rpri,R,M,N,Ext,x);
        
        DL_glc(m) = abs(alpha0/belta0);
        [SINR_GLC(m),SOI_GLC(m)] = PerformanceF(EnerS,Q,Rglc,R,1,ad,adS);
        
        % SVE of GLC
        ad_glc_sve = SVEfunction(Ang_left,Ang_right,Ang_interval,kd, Rglc);
        [SINR_GLC_SVE(m),SOI_GLC_SVE(m)] = PerformanceF(EnerS,Q,Rglc,R,1,ad_glc_sve,adS);
        
        %% method of HKB             
%         ad = (exp(1j * kd * sin(DOAS + delta) * (1 : M)'));
        [rhoHKB, Rhkb] = HKBfunction(ad, R, M);
        DL_hkb(m) = abs(rhoHKB);
        [SINR_HKB(m),SOI_HKB(m)] = PerformanceF(EnerS,Q,Rhkb,R,1,ad,adS);
        
        % applying the SVE into HKB
%         ad_est = SVEfunction(Ang_left,Ang_right,Ang_interval,kd, R);
%         [rhoHKB_SVE, Rhkb_SVE] = HKBfunction(ad_est, R, M);
        [SINR_HKB_SVE(m),SOI_HKB_SVE(m)] = PerformanceF(EnerS,Q,Rhkb,R,1,ad_est,adS);
        
        %% method of DL-CGLS        
        ad_norm = ad/sqrt(M);
        DLcg = ad_norm' * R * ad_norm; % the empirical DL of CG method
              
        % LCMV beamforming
        ConsMat = ad;
        Gv = 1;
        SCtype = 'GDP'; %Stopping criteria type: GDP for generalized discrepancy principle; EE for error estimation; RVE for Ritz value estimation
        tol_dB = 5;
        [Wcg] = CGLSfunction(ConsMat, Gv, x, tol_dB, DLcg, SCtype);
       
        DL_cg(m) = DLcg;
        SINR_cg(m) = 10 * log10(10^(EnerS/10) * abs(Wcg' * adS)^2 / real(Wcg' * Q * Wcg));
        SOI_cg(m) = 10 * log10(abs(Wcg' * R0 * Wcg));
        
        % considering the SVE
        ConsMat = ad_est;
        [Wcg] = CGLSfunction(ConsMat, Gv, x, tol_dB, DLcg, SCtype);
        SINR_cg_sve(m) = 10 * log10(10^(EnerS/10) * abs(Wcg' * adS)^2 / real(Wcg' * Q * Wcg));
        SOI_cg_sve(m) = 10 * log10(abs(Wcg' * R0 * Wcg));
        
        
        %% method of QC
        Ang_left_qc = DOAS - 2*pi/180;
        Ang_right_qc = DOAS + 2*pi/180;
        TestPoint = 4;     
        DLqc = 1;
        if sin(Ang_right_qc) - sin(Ang_left_qc) <= 2/M
            [Wqc, DLqc] = QuadConsFunction(Ang_left_qc, Ang_right_qc, TestPoint, kd, M, R, DLqc);
        else
            disp('warnning: The angle threshold of QC method is too large')
            break;
        end

        DL_qc(m) = DLqc;
        SINR_qc(m) = 10 * log10(10^(EnerS/10) * abs(Wqc' * adS)^2 / real(Wqc' * Q * Wqc));
        SOI_qc(m) = 10 * log10(abs(Wqc' * R0 * Wqc));
        
        
    end
     % DL
    DL_glc_ls_list(NIdx) = mean(DL_glc_ls);
    DL_glc_mwf_list(NIdx) = mean(DL_glc_mwf);
    DL_glc_list(NIdx) = mean(DL_glc);
    DL_hkb_list(NIdx) = mean(DL_hkb);
    DL_cg_list(NIdx) = mean(DL_cg);
    DL_qc_list(NIdx) = mean(DL_qc);
    
    % SINR
    SINR_opt_list(NIdx) = mean(SINR_opt);
    SINR_GLC_LS_list(NIdx) = mean(SINR_GLC_LS);
    SINR_GLC_LS_SVE_list(NIdx) = mean(SINR_GLC_LS_SVE);    
    SINR_GLC_MWF_list(NIdx) = mean(SINR_GLC_MWF);
    SINR_GLC_MWF_SVE_list(NIdx) = mean(SINR_GLC_MWF_SVE);  
    SINR_GLC_list(NIdx) = mean(SINR_GLC);
    SINR_GLC_SVE_list(NIdx) = mean(SINR_GLC_SVE);    
    SINR_HKB_list(NIdx) = mean(SINR_HKB);
    SINR_HKB_SVE_list(NIdx) = mean(SINR_HKB_SVE);    
    SINR_CG_list(NIdx) = mean(SINR_cg);
    SINR_CG_SVE_list(NIdx) = mean(SINR_cg_sve);   
    SINR_QC_list(NIdx) = mean(SINR_qc);
    SINR_SCB_list(NIdx) = mean(SINR_SCB);
    
    % SOI
    SOI_opt_list(NIdx) = mean(SOI_opt);    
    SOI_GLC_LS_list(NIdx) = mean(SOI_GLC_LS);
    SOI_GLC_LS_SVE_list(NIdx) = mean(SOI_GLC_LS_SVE);    
    SOI_GLC_MWF_list(NIdx) = mean(SOI_GLC_MWF);
    SOI_GLC_MWF_SVE_list(NIdx) = mean(SOI_GLC_MWF_SVE);   
    SOI_GLC_list(NIdx) = mean(SOI_GLC);
    SOI_GLC_SVE_list(NIdx) = mean(SOI_GLC_SVE);    
    SOI_HKB_list(NIdx) = mean(SOI_HKB);
    SOI_HKB_SVE_list(NIdx) = mean(SOI_HKB_SVE);    
    SOI_CG_list(NIdx) = mean(SOI_cg);
    SOI_CG_SVE_list(NIdx) = mean(SOI_cg_sve);    
    SOI_QC_list(NIdx) = mean(SOI_qc);
    SOI_SCB_list(NIdx) = mean(SOI_SCB);
end

%% figures of the DL value
figure(1),
plot(N_list, abs(DL_glc_list),'r','linewidth',2)
set(gca,'FontName','Times New Roman','FontSize',20);
xlabel('Number of Snapshots','FontName','Times New Roman','FontSize',20);
ylabel({'{\rho}_{glc}'},'FontName','Times New Roman','FontSize',20);

figure(2),
plot(N_list, abs(DL_hkb_list),'r','linewidth',2)
set(gca,'FontName','Times New Roman','FontSize',20);
xlabel('Number of Snapshots','FontName','Times New Roman','FontSize',20);
ylabel({'{\rho}_{hkb}'},'FontName','Times New Roman','FontSize',20);

figure(3),
plot(N_list, abs(DL_cg_list),'r','linewidth',2)
set(gca,'FontName','Times New Roman','FontSize',20);
xlabel('Number of Snapshots','FontName','Times New Roman','FontSize',20);
ylabel({'{\rho}_{cg}'},'FontName','Times New Roman','FontSize',20);

figure(4)
plot(N_list, abs(DL_qc_list), 'r','linewidth',2)
set(gca,'FontName','Times New Roman','FontSize',20);
xlabel('Number of Snapshots','FontName','Times New Roman','FontSize',20);
ylabel({'{\rho}_{qc}'},'FontName','Times New Roman','FontSize',20);

figure(5),
plot(N_list,abs(DL_glc_ls_list),'r','linewidth',2)
set(gca,'FontName','Times New Roman','FontSize',20);
xlabel('Number of Snapshots','FontName','Times New Roman','FontSize',20);
ylabel({'{\rho}_{ls}'},'FontName','Times New Roman','FontSize',20);

figure(6),
plot(N_list,abs(DL_glc_mwf_list),'r','linewidth',2)
set(gca,'FontName','Times New Roman','FontSize',20);
xlabel('Number of Snapshots','FontName','Times New Roman','FontSize',20);
ylabel({'{\rho}_{mwf}'},'FontName','Times New Roman','FontSize',20);

%% figures of the N-SINR curves with SVE modification
RGB_table = [ {[0 0 0]}; {[0 255 255]}; {[0 205 205]};
    {[84 255 159]}; {[0 205 0]}; {[132 112 255]}; 
    {[0 0 205]}; {[105 105 105]}; {[255 0 0]};
    {[255 99 71]}; {[255 0 205]}; {[148 0 211]}; {[0 0 255]}];

figure(70),
h = plot(N_list,SINR_opt_list,'k--',...
    N_list, SINR_GLC_list,'->',...
    N_list, SINR_GLC_SVE_list,'-.>',...
    N_list,SINR_HKB_list,'-o',...
    N_list,SINR_HKB_SVE_list,'-.o',...
    N_list, SINR_CG_list, '-v',...
    N_list, SINR_CG_SVE_list, '-.v',...
    N_list, SINR_QC_list,'-p',...
    N_list,SINR_GLC_LS_list,'-s',...
    N_list,SINR_GLC_LS_SVE_list,'-.s',...
    N_list,SINR_GLC_MWF_list,'-<',...
    N_list,SINR_GLC_MWF_SVE_list,'-.<',...
    N_list, SINR_SCB_list,'-d');
set(h, 'linewidth',2.5,'MarkerSize',6)
for idx = 1 : length(h)
    set(h(idx),'color',cell2mat(RGB_table(idx))/255);
end
xlabel('Number of Snapshots','FontName','Times New Roman','FontSize',20);
ylabel({'{SINR\out}','dB'},'FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20); 
hleg = legend('SINR_{opt}','GLC','GLC-SVE','HKB','HKB-SVE','CG','CG-SVE','QC','LS-TMMSE','LS-TMMSE-SVE','MWF-TMMSE','MWF-TMMSE-SVE','SCB','Location','SouthEast');
set(hleg,'Fontsize',10);
axis([N_list(1),N_list(end),min(SINR_SCB_list)-5,max(SINR_opt_list)+2])
set (gcf,'Position',[10,10,1200,900])
nummarkers(h, 10);

%% figures of the N-SINR curves without SVE modification
nplots = [1 2 4 6 8 9 11 13];
color = cell(length(nplots), 1);
for idx = 1 : length(nplots)
    color{idx} = get(h(nplots(idx)),'color');
end

figure(80),
h_elimsve = plot(N_list,SINR_opt_list,'k--',...
    N_list, SINR_GLC_list,'->',...
    N_list,SINR_HKB_list,'-o',...
    N_list, SINR_CG_list, '-v',...
    N_list, SINR_QC_list,'-p',...
    N_list,SINR_GLC_LS_list,'-s',...
    N_list,SINR_GLC_MWF_list,'-<',...
    N_list, SINR_SCB_list,'-d');
set(h_elimsve, 'linewidth',2.5,'MarkerSize',6)
lh_elimsve = length(h_elimsve);
for n = 1 : lh_elimsve
    set(h_elimsve(n), 'color',color{n});
end
xlabel('Number of Snapshots','FontName','Times New Roman','FontSize',20);
ylabel({'{SINR\out}','dB'},'FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20); 
hleg = legend('SINR_{opt}','GLC','HKB','CG','QC','LS-TMMSE','MWF-TMMSE','SCB','Location','SouthEast');
set(hleg,'Fontsize',10);
axis([N_list(1),N_list(end),min(SINR_SCB_list)-5,max(SINR_opt_list)+2])
set (gcf,'Position',[10,10,1200,900])
nummarkers(h_elimsve, 10);
