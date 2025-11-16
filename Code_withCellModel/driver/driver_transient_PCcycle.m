% driver_run_eqPCcycle.m
%   Driver script to run biogeochemical model, without optimizing parameters.
%   This run only solves the phosphorus and carbon cycles, using 
%   predetermined parameters
%
%   for Sullivan 2025 paper: 
%   AIM: to quantify the impact of productivity hangover. 
%   Uses the reoptimized Nature2023 model parameters from 
%   reoptNature_with_dop_GM15_npp1_CTL_He_PCO_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00.mat
%
%   Increase surface production to simulate idealized nutrient fertilization
%   experiments, and track the transient response of the ocean carbon cycle.
% ------------------------------------------------------------------------
%clc; clear all; close all
global iter
iter = 0 ;
on   = true  ; off  = false ;
format short

% --- addpath to model code -----
addpath('../src/')
%addpath('../src_reoptNature/')

% test1_eqPcycle_with_DOPl_gamma1pct_from_reoptNature_with_dop_GM15_npp1

VerName = 'transient_from_reoptNature_with_dop_GM15_npp1_'; 		% optional version name. leave as an empty character array
					% or add a name ending with an underscore
VerNum = '';		% optional version number for testing

% Choose C2P function
par.C2Pfunctiontype = 'P';
% 'P' -> PO4 function ; 'C' -> Cell model; 'T' -> Temperature function; 'R' -> constant value (Redfield)
% 
par.nppVer = 1; % 1 = CbPM; 2 = CAFE; (Nowicki)
%
GridVer  = 91  ;
operator = 'A' ;
% GridVer: choose from 90 and 91; Ver 90 is for a Transport
% operator without diapycnal mixing but optimized using DIP ;
% Ver 91 include a bunch of operators that include diapycnal
% mixing. These operators represent sensiviity tests on He
% constraint and on mixing parameterizations (DeVries et al, 2018).
% A -> CTL_He; B -> CTL_noHe; C -> KiHIGH_He; D -> KiHIGH_noHe;
% E -> KvHIGH_KiLOW_He; F -> KvHIGH_KiLOW_noHe; G -> KiLOW_He;
% H -> KiLOW_noHe; I -> KvHIGH_He; J -> KvHIGH_noHe; K -> KvHIGH_KiHIGH_noHe
% L -> CTL_He_48layer

par.nl = 2; % number of layers in the model euphotic zone (doesn't change)

Gtest = off ; 
Htest = off ;
par.optim   = off ; 
par.Cmodel  = on ; 
par.Omodel  = off ; 
par.Simodel = off ;
par.Cisotope  = off  ;
par.dynamicP = off ; % if on, cell model uses modeled DIP. if off, cell model uses WOA observed DIP field.
par.LoadOpt = on ; % if load optimial parameters. 
% to load parameter values from a run with a different name.
%par.fxhatload = '../../output/optPonly_CTL_He_P_xhat.mat';
%par.fxhatload = '/DFS-L/DATA/primeau/hojons1/Nature2023_BGC_reoptimized/src_Nature_parameter_Megan/MSK91/CTL_He_PCO_Gamma0_kl12h_O5_POC2DIC_GM15_Nowicki_npp1_aveTeu_diffSig_O2C_uniEta_DICrmAnthro_2L_Pnormal_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00_xhat.mat' 
%par.fxhatload = '../output/test_equalbkCP_reoptNature_with_dop_GM15_npp1_CTL_He_xhat.mat'
par.fxhatload = '../output/PNAS2025/reoptNature_with_dop_GM15_npp1_CTL_He_PCO_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00_xhat.mat'

% to use different model output for initial CO guess. 
%par.fnameload = '/DFS-L/DATA/primeau/hojons1/Nature2023_BGC_reoptimized/src_Nature_parameter_Megan/MSK91/CTL_He_PCO_Gamma0_kl12h_O5_POC2DIC_GM15_Nowicki_npp1_aveTeu_diffSig_O2C_uniEta_DICrmAnthro_2L_Pnormal_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00.mat' ;
par.fnameload = '../output/PNAS2025/reoptNature_with_dop_GM15_npp1_CTL_He_PCO_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00.mat'

% to load presaved input data fields in par
% par.fparload = '../output/PNAS2025/reoptNature_with_dop_GM15_npp1_CTL_He_PCO_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00_par.mat'

par.dopscale = 1.0 ;
par.dipscale = 1.0 ;
par.dicscale = 1.0 ;
par.docscale = 1.0 ; % factor to weigh DOC in the objective function
par.alkscale = 1.0 ;
par.o2scale  = 1.0 ;
% P model parameters
par.opt_sigP  = on ; 
par.opt_Q10P  = on ;
par.opt_kdP   = on ;
par.opt_bP_T  = on ; 
par.opt_bP    = on ;
par.opt_alpha = on ;
par.opt_beta  = on ;
% C model parameter
par.opt_sigC  = on ; 
par.opt_kru   = on ;
par.opt_krd   = on ;
par.opt_etau  = on ;
par.opt_etad  = off ; %keep off
par.opt_bC_T  = on ;
par.opt_bC    = on ; 
par.opt_d     = on ;
par.opt_Q10C  = on ;
par.opt_kdC   = on ; 
par.opt_R_Si  = on ; 
par.opt_rR    = on ; 
% --- C:P function parameters -----
% phosphate-dependent function parameters
par.opt_cc    = off ;
par.opt_dd    = off ; 
% temperature-dependent function parameters
par.opt_ccT   = on ; 
par.opt_ddT   = on ;
% Trait-based Cellular Growth Model parameters
par.opt_Q10Photo     = on ; % opt
par.opt_fStorage     = on ; % opt
par.opt_fRibE 	     = off ; 
par.opt_kST0 	     = on ; % opt
par.opt_PLip_PCutoff = off ;
par.opt_PLip_scale   = off ;
par.opt_PStor_rCutoff = on; % opt
par.opt_PStor_scale  = off ;
par.opt_alphaS       = on ; % opt
par.opt_gammaDNA	 = off ;
% O model parameters
par.opt_O2C_T = off ;
par.opt_rO2C  = on ;
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = off ;
par.opt_bt    = on  ;
par.opt_aa    = on  ;
par.opt_bb    = on  ;
%
%-------------load data and set up parameters---------------------
SetUp ;                      

% save results 
% ATTENTION: Change this directory to where you want to
% save your output files
output_dir = sprintf('../output/PNAS2025_transient/'); 

if ~isdir(output_dir)
    command = strcat('mkdir', " ", output_dir) ;
    system(command) ;
end

VER = strcat(output_dir,VerName,TRdivVer);
catDOC = ''; % sprintf('_DOC%0.2g_DOP%0.2g',par.docscale,par.dopscale); % used to add scale factors to file names
% Create output file names based on which model(s) is(are) optimized
%if Gtest == on
%    fname = strcat(VER,'_GHtest');
%elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off & par.Cellmodel == off)
        fname = strcat(VER,'_P',VerNum);
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off & par.Cellmodel == off)
        base_name = strcat(VER,'_PC',VerNum);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off & par.Cellmodel == off)
        base_name = strcat(VER,'_PCO',VerNum);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on & par.Cellmodel == off)
        base_name = strcat(VER,'_PCSi',VerNum);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on & par.Cellmodel == off)
        base_name = strcat(VER,'_PCOSi',VerNum);
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == off & par.Omodel == off & par.Simodel == off & par.Cellmodel == on) % cell model does nothing if C model is not on, so this case =Ponly
        base_name = strcat(VER,'_PCell',VerNum);
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off & par.Cellmodel == on)
        base_name = strcat(VER,'_PCCell',VerNum);
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off & par.Cellmodel == on)
		base_name = strcat(VER,'_PCOCell',VerNum);
		fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on & par.Cellmodel == on)
        base_name = strcat(VER,'_PCOSiCell',VerNum);
        fname = strcat(base_name,catDOC);
    end
%end

% -------------------- Set up output files ---------------
% -------------------- Set up output files ---------------
par.fname = strcat(fname,'.mat') ;
fxhat     = strcat(fname,'_xhat.mat') 
fxpar     = strcat(fname,'_par.mat');
par.fxhat = fxhat ;
if Htest ==on
	fGHtest = strcat(fname,'_GHtest.mat')  ;
end
par.fxhat = fxhat ;
par.fxpar = fxpar ;

% -------------------update initial guesses --------------
if isfile(par.fname)
    fprintf('loading initial guess on C and O from file: %s \n',par.fname)
    load(par.fname,'data')
end 

% -------------------update initial guesses --------------
if isfile(par.fnameload)
    fprintf('loading initial guess on C and O from file: %s \n',par.fnameload)
    load(par.fnameload,'data')
end 

% to timestep P, need to define initial guess as optimal model solution

%---------------- inital guesses on C and O ---------------
if par.Cmodel == on 
    pco2atm = par.pco2_air(1) ;  % uatm
    GC  = [data.DIC(iwet); data.POC(iwet); data.DOC(iwet); data.PIC(iwet); ...
           data.ALK(iwet); data.DOCl(iwet); data.DOCr(iwet); pco2atm];
    % GC  = real(GC) + 1e-6*randn(7*nwet+1,1) ;

    % Set total carbon in System (Atm + Ocn)
	%par.TotalCset = 6.2816e+18 ; %[mol C]
	Natm    = 1.773e20  ; %  molar volume of atmosphere
	par.Natm = Natm ;
	TC_ocean = ((data.DIC(iwet)+data.DOC(iwet)+data.POC(iwet)+data.PIC(iwet)+data.DOCl(iwet)+data.DOCr(iwet))'*par.dVt(iwet)).*1e-3; %[molC]
	par.totalCarbon = (data.DIC(iwet)+data.DOC(iwet)+data.POC(iwet)+data.PIC(iwet)+data.DOCl(iwet)+data.DOCr(iwet))'*par.dVt(iwet).*1e-3 + pco2atm*Natm*1e-6  ;  %[molC]
	fprintf('Initial pco2atm = %4.2f ppm ; totalC_ocean = %6.5e mol C; totalCarbon = %6.5e mol C \n', pco2atm, TC_ocean, par.totalCarbon) ;
	clear TC_ocean

end 
if par.Omodel == on 
    GO  = real(data.O2(iwet)) ; % + 1e-9*randn(par.nwet,1);
end 

% save initial solution in par?
par.DIP = data.DIP(iwet); 
par.POP = data.POP(iwet);
par.DOP = data.DOP(iwet);
par.DOPl= data.DOPl(iwet);
par.DIC = data.DIC(iwet);
par.POC = data.POC(iwet);
par.DOC = data.DOC(iwet);
par.PIC = data.PIC(iwet);
par.ALK = data.ALK(iwet);
par.DOCl= data.DOCl(iwet);
par.DOCr= data.DOCr(iwet);
par.O2  = data.O2(iwet);

%% 
%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPar(par)  ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

x0    = par.p0 ;
nip = length(x0);
%-------------------set up timestepper -------------------------
spa = par.spa;


% Default schedule
%                      1 hr    1 day   1 month  1 year   4 years   
dt_size =        spa./[365*24  365     12       1        0.25   ];  % step sizes in seconds
nsteps =              [ 50      50      50       20       30     ]; % number of steps for each size
par.dt_schedule = dt_size ;
par.nsteps_schedule = nsteps ;

total_nsteps = sum(nsteps);
total_time_yr = sum(dt_size .* nsteps) / spa;
fprintf('Time-step Schedule: %d total steps over %.2f years\n', total_nsteps, total_time_yr);


% allocate output  
% diagnostics to save at every timestep: 
% pco2atm, volume integrated DIC, integrated PNPP, CNPP
Tout = zeros(1,total_nstep);
pco2atmout = zeros(1,total_nstep);

t0 = 0; 
  t = t0;
  Tout(1) = t;
  fprintf('Before time stepping method (Steady-state solutions in pre-industrial era)...\n');
  fprintf('...Time: %4.2f, AtmCO2: %4.2f uatm,  avgDIC: %7.6g mmol/m3\n', ...
    t,par.pco2atm,sum(par.DIC.*par.dVt(iwet))/sum(par.dVt(iwet))); % mean(par.DIC)

    totalDIC = sum(par.DIC.*par.dVt(iwet).*1e-3); % units = mol C
    totalCO2atm = par.pco2atm * par.Natm * 1e-6; % mol C
    fprintf('...Ocn: avgDIC = %7.4f mmol/m3 \n',sum(par.DIC.*par.dVt(iwet))/sum(par.dVt(iwet)));
    fprintf('...Ocn: Integrated total DIC = %10.3e Pg C \n',totalDIC*12*1e-15)
    fprintf('...Atm: Integrated total CO2 = %10.3e Pg C \n',totalCO2atm*12*1e-15);



% 3) Time integration with scheduled stepping
% -----------------------
% trapezoid + euler forward
    % dX/dt + F(X,t) =0 where F(X,t) = J*X + f(X,t)
    % 
    % X(n+1) - X(n) + (dt/2)*J*(X(n+1)+X(n)) + dt*fn = 0
    % ==> A*X(n+1) = B*X(n) - dt*f(n)
    % X(n+1) = A\(B*X(n)-dt*f(n))

% t = t0;
% global_step = 1; % global step counter for saving
% for dt_idx = 1:length(dt_size)
%     dt = dt_size(dt_idx);
%     n_substeps = nsteps(dt_idx);

%     for ii = 1:n_substeps 
%         global_step = global_step + 1; 
%     end
% end


%% Time step model
% Set up
% solve/load steady state solution 
% check the data and par for steady-state drive the time stepping method 
    Xin.P = [par.DIP;par.POP;par.DOP;par.DOPl];
    Xin.C = [par.DIC;par.POC;par.DOC;par.PIC;par.ALK;par.DOCl;par.DOCr; par.pco2atm];
    Xin.O2 = [par.O2];
% set up iteration counter 
t = t0;
global_step = 1; % global step counter for saving
for dt_idx = 1:length(dt_size)
    dt = dt_size(dt_idx);
    n_substeps = nsteps(dt_idx);
    % load previous iteration solution
    Xin.P = [par.DIP;par.POP;par.DOP;par.DOPl];
    Xin.C = [par.DIC;par.POC;par.DOC;par.PIC;par.ALK;par.DOCl;par.DOCr; par.pco2atm];
    Xin.O2 = [par.O2];

    % set up time stepper for P and C (only need to factor trapezoid matrix once for each step size dt)
    % run Peqn
    [F_P,J_P,par] = Peqn(Xin.P, par); 
    % Evaluate RHS and Jacobian at current state (X, t)
    % Build trapezoidal matrices
    I_P  = speye(numel(Xin.P(:)));
    A_P  = I_P + 0.5*dt*J_P;
    B_P  = I_P - 0.5*dt*J_P;
    A_Pfactored = mfactor(A_P); 

    % run CeqnAtm
    [F_C,J_C,par] = CeqnAtm(Xin.C, par);
    % Evaluate RHS and Jacobian at current state (X, t)
    % Build trapezoidal matrices
    I_C  = speye(numel(Xin.C(:)));
    A_C  = I_C + 0.5*dt*J_C;
    B_C  = I_C - 0.5*dt*J_C;
    A_Cfactored = mfactor(A_C);

    % start iterating through n_substeps for each dt
    for ii = 1:n_substeps 
        global_step = global_step + 1; 

        % load previous iteration solution
        Xin.P = [par.DIP;par.POP;par.DOP;par.DOPl];
        Xin.C = [par.DIC;par.POC;par.DOC;par.PIC;par.ALK;par.DOCl;par.DOCr; par.pco2atm];
        Xin.O2 = [par.O2];

        % run Peqn
        [F_P,J_P,par] = Peqn(Xin.P, par); 
        % Right-hand side
        rhs_P = B_P*Xin.P - dt*F_P;
        Xout.P  = mfactor(A_Pfactored, rhs_P);        

        % update par with Peqn solution
        par.DIP = Xout.P(1:nwet);
        par.POP = Xout.P(1*nwet+1:2*nwet);
        par.DOP = Xout.P(2*nwet+1:3*nwet);
        par.DOPl= Xout.P(3*nwet+1:4*nwet);

        % run CeqnAtm
        [F_C,J_C,par] = CeqnAtm(Xin.C, par);
        % Right-hand side
        rhs_C = B_C*Xin.C - dt*F_C;
        Xout.C  = mfactor(A_Cfactored, rhs_C);
        % update par with CeqnAtm solution
        par.DIC    = Xout.C(1:nwet);
        par.POC    = Xout.C(1*nwet+1:2*nwet);
        par.DOC    = Xout.C(2*nwet+1:3*nwet);
        par.PIC    = Xout.C(3*nwet+1:4*nwet);
        par.ALK    = Xout.C(4*nwet+1:5*nwet);
        par.DOCl   = Xout.C(5*nwet+1:6*nwet);
        par.DOCr   = Xout.C(6*nwet+1:7*nwet);
        par.pco2atm= Xout.C(7*nwet+1);  

        % save output at each time step


    end
end



% %%%%%%%%%%%%%%%%%%   steady state BGC model   %%%%%%%%%%%%%%%%%%%%%%%%
function [f,J] =PCeqns(x0, par)
    clear data
    x = x0;
    iter = 0;
    %[f,fx,fxx,data] = neglogpost(xsol,par);
    %fprintf('----neglogpost complete----\n')
    fprintf('\ncurrent time is:      %s\n',datetime('now')) ;
    fprintf('current iteration is: %d \n',iter) ;

    % print and save current parameter values to
    % a file that is used to reset parameters ;
    PrintPar(x, par) ;    
    % increment iteration counter
    iter = iter + 1  ;

    nx   = length(x) ; % number of parameters
    dVt  = par.dVt   ;
    M3d  = par.M3d   ;
    iwet = par.iwet  ;
    nwet = par.nwet  ;
    %
    f    = 0 ;
    %%%%%%%%%%%%%%%%%%   Solve P    %%%%%%%%%%%%%%%%%%%%%%%%
    idip = find(par.po4raw(iwet) > 0.05) ;
    Wp   = d0(dVt(iwet(idip))/sum(dVt(iwet(idip)))) ;
    mu   = sum(Wp*par.po4raw(iwet(idip)))/sum(diag(Wp)) ;
    var  = sum(Wp*(par.po4raw(iwet(idip))-mu).^2)/sum(diag(Wp)) ;
    Wip  = par.dipscale*Wp/var ;

    idop = find(par.dopraw(iwet) > 0.0) ;
    Wp   = d0(dVt(iwet(idop))/sum(dVt(iwet(idop)))) ;
    mu   = sum(Wp*par.dopraw(iwet(idop)))/sum(diag(Wp)) ;
    var  = sum(Wp*(par.dopraw(iwet(idop))-mu).^2)/sum(diag(Wp)) ;
    Wop  = par.dopscale*Wp/var ;
    %
    %tic 
    [par, P, Px, Pxx] = eqPcycle(x, par) ;
    DIP  = M3d+nan  ;  DIP(iwet)  = P(1+0*nwet:1*nwet) ;
    POP  = M3d+nan  ;  POP(iwet)  = P(1+1*nwet:2*nwet) ;
    DOP  = M3d+nan  ;  DOP(iwet)  = P(1+2*nwet:3*nwet) ;
    DOPl = M3d+nan  ;  DOPl(iwet) = P(1+3*nwet:4*nwet) ;
    %toc 
    par.Px   = Px  ;
    par.Pxx  = Pxx ;
    par.DIP  = DIP(iwet) ; %need to update DIP field for C cycle run
    par.POP  = POP(iwet) ;
    par.DOP  = DOP(iwet) ;
    par.DOPl = DOPl(iwet) ;

    data.DIP = DIP ; data.POP  = POP  ;
    data.DOP = DOP ; data.DOPl = DOPl ;
    % DIP & DOP error
    DOP = DOP + DOPl; % sum of semilabile and labile DOP ;
    eip = DIP(iwet(idip)) - par.po4raw(iwet(idip)) ;
    eop = DOP(iwet(idop)) - par.dopraw(iwet(idop)) ;
    f  = f + 0.5*(eip.'*Wip*eip) + 0.5*(eop.'*Wop*eop); 
    f_components.DIP = 0.5*(eip.'*Wip*eip);
    f_components.DOP = 0.5*(eop.'*Wop*eop); 

    
    
    %%%%%%%%%%%%%%%%%%   End Solve P    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%     Solve C   %%%%%%%%%%%%%%%%%%%%%%%%
    if (par.Cmodel == on)
        idic = find(par.dicraw(iwet) > 0) ;
        Wic  = d0(dVt(iwet(idic))/sum(dVt(iwet(idic)))) ;
        mu   = sum(Wic*par.dicraw(iwet(idic)))/sum(diag(Wic)) ;
        var  = sum(Wic*(par.dicraw(iwet(idic))-mu).^2)/sum(diag(Wic));
        Wic  = par.dicscale*Wic/var  ;
        
        ialk = find(par.alkraw(iwet) > 0) ;
        Wlk  = d0(dVt(iwet(ialk))/sum(dVt(iwet(ialk)))) ;
        mu   = sum(Wlk*par.alkraw(iwet(ialk)))/sum(diag(Wlk)) ;
        var  = sum(Wlk*(par.alkraw(iwet(ialk))-mu).^2)/sum(diag(Wlk));
        Wlk  = par.alkscale*Wlk/var  ;
        
        idoc = find(par.docraw(iwet) > 0) ;
        Woc  = d0(dVt(iwet(idoc))/sum(dVt(iwet(idoc)))) ;
        mu   = sum(Woc*par.docraw(iwet(idoc)))/sum(diag(Woc)) ;
        var  = sum(Woc*(par.docraw(iwet(idoc))-mu).^2)/sum(diag(Woc));
        Woc  = par.docscale*Woc/var ;
        %tic 
        [par, C, Cx, Cxx] = eqCcycleAtm(x, par) ;
        DIC  = M3d+nan ;  DIC(iwet)  = C(0*nwet+1:1*nwet) ;
        POC  = M3d+nan ;  POC(iwet)  = C(1*nwet+1:2*nwet) ;
        DOC  = M3d+nan ;  DOC(iwet)  = C(2*nwet+1:3*nwet) ;
        PIC  = M3d+nan ;  PIC(iwet)  = C(3*nwet+1:4*nwet) ;
        ALK  = M3d+nan ;  ALK(iwet)  = C(4*nwet+1:5*nwet) ;
        DOCl = M3d+nan ;  DOCl(iwet) = C(5*nwet+1:6*nwet) ;
        DOCr = M3d+nan ;  DOCr(iwet) = C(6*nwet+1:7*nwet) ;
        pco2atm = C(7*nwet+1);
       % toc

        par.DIC = DIC(iwet) ;
        par.POC = POC(iwet) ;
        par.DOC = DOC(iwet) ;
        par.DOCl = DOCl(iwet) ;
        par.DOCr = DOCr(iwet) ;
        par.pco2atm = pco2atm ;
        % DIC = DIC + par.dicant  ;
        par.Cx    = Cx   ;  par.Cxx   = Cxx ;
        data.DIC  = DIC  ;  data.POC  = POC ;
        data.DOC  = DOC  ;  data.PIC  = PIC ;
        data.ALK  = ALK  ;  data.DOCr = DOCr ;
        data.DOCl = DOCl ;  data.pco2atm = pco2atm ;
        try
            data.C2P = M3d+nan ; data.C2P(iwet) = par.C2P; 
        catch ME
            fprintf('error in %s (line %d): %s \n', ME.stack(1).name,ME.stack(1).line,ME.message);
            fprintf('Unable to store C2P in data struct. \n');
        end
        % DOC error
        DOC = DOC + DOCr + DOCl; % sum of labile and refractory DOC ;
        eic = DIC(iwet(idic)) - par.dicraw(iwet(idic)) ;
        eoc = DOC(iwet(idoc)) - par.docraw(iwet(idoc)) ;
        elk = ALK(iwet(ialk)) - par.alkraw(iwet(ialk)) ;
        f   = f + 0.5*(eic.'*Wic*eic) + 0.5*(eoc.'*Woc*eoc) + ...
              0.5*(elk.'*Wlk*elk);

        f_components.DIC = 0.5*(eic.'*Wic*eic);
        f_components.DOC = 0.5*(eoc.'*Woc*eoc);
        f_components.ALK = 0.5*(elk.'*Wlk*elk);

        % % print carbon system info
	    % fprintf('Atm CO2 concentration: %3.2f ppm \n', pco2atm);
		% %DIC
		% DICtmp = model.DIC.*par.dVt.*1e-3;  %units = mol C
    	% tDIC = sum(DICtmp(iwet),'all');
    	% fprintf('Integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC,tDIC*12*1e-15)
        % % all forms of carbon
        % DICtmp = (model.DIC+model.DOC + model.POC +model.PIC).*par.dVt.*1e-3;
        % tC_ocn = sum(DICtmp(iwet),'all');
		% tC_all = tC_ocn + model.pco2atm*par.Natm*1e-6;
		% fprintf('Atm pCO2 from prescribed TC - TC_ocean: %3.2f ppm \n',(par.totalCarbon - tC_ocn)/par.Natm*1e6) ;
		% fprintf('Integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n',tC_ocn,tC_ocn*12*1e-15);
		% fprintf('total C in System (Atm+Ocn) =   %10.4e mol C   (%6.3e Pg C) \n\n',tC_all,tC_all*12*1e-15);
		% model.totalC = tC_all;
		% model.totalC_ocean = tC_ocn;
    end
    %%%%%%%%%%%%%%%%%%   End Solve C    %%%%%%%%%%%%%%%%%%%

    % Print and save objective function subcomponent values
    data.f = f;
    data.f_components = f_components;

    fprintf('current objective function value is: %3.3e \n\n',f) 
    fprintf('current objective function value for fit to DIP is %3.3e \n',f_components.DIP) 
    fprintf('current objective function value for fit to DOP is %3.3e \n',f_components.DOP) 
    if (par.Cmodel == on)
        fprintf('current objective function value for fit to DIC is %3.3e \n',f_components.DIC) 
        fprintf('current objective function value for fit to DOC is %3.3e \n',f_components.DOC) 
        fprintf('current objective function value for fit to ALK is %3.3e \n',f_components.ALK) 
    end


    %% note: skipping save for testing
    if exist(par.fname, 'file')
        reply = input(sprintf('WARNING: File ( %s ) already exists. \nDo you want to overwrite this file? Y/N: ', par.fname), 's');
        if strcmpi(reply, 'Y')
            fprintf('Overwriting File... \n');
            fprintf('saving model solution to file: %s \n',par.fname)
            save(par.fname, 'data')
        else
            fprintf('Execution stopped by User.\n');
            fprintf('--------------------------\n\n');
            return;
        end
    else
       fprintf('saving model solution to file: %s \n',par.fname)
       save(par.fname, 'data')
    end


fprintf('-------------- end! ---------------\n');
end