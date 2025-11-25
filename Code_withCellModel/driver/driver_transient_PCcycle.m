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

VerName = 'transient_test_OIF_PC_noAtm_from_reoptNature_with_dop_GM15_npp1_'; 		% optional version name. leave as an empty character array
					% or add a name ending with an underscore
VerNum = '';		% optional version number for testing

par.VerName = VerName;

par.equalbkCP_flag = false ; % reset parameters to make remin rates equal for P and C 
par.fertilize_flag = true; % if true, increase surface productivity for 1 year

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
par.saveall = false; 
par.save_stride = 12;   % monthly saves by default

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
par.opt_sigP  = off ; 
par.opt_Q10P  = off ;
par.opt_kdP   = off ;
par.opt_bP_T  = off ; 
par.opt_bP    = off ;
par.opt_alpha = off ;
par.opt_beta  = on ;
% C model parameter
par.opt_sigC  = off ; 
par.opt_kru   = off ;
par.opt_krd   = off ;
par.opt_etau  = off ;
par.opt_etad  = off ; %keep off
par.opt_bC_T  = off ;
par.opt_bC    = off ; 
par.opt_d     = off ;
par.opt_Q10C  = off ;
par.opt_kdC   = off ; 
par.opt_R_Si  = off ; 
par.opt_rR    = off ; 
% --- C:P function parameters -----
% phosphate-dependent function parameters
par.opt_cc    = off ;
par.opt_dd    = off ; 
% temperature-dependent function parameters
par.opt_ccT   = off ; 
par.opt_ddT   = off ;
% Trait-based Cellular Growth Model parameters
par.opt_Q10Photo     = off ; % opt
par.opt_fStorage     = off ; % opt
par.opt_fRibE 	     = off ; 
par.opt_kST0 	     = off ; % opt
par.opt_PLip_PCutoff = off ;
par.opt_PLip_scale   = off ;
par.opt_PStor_rCutoff = off; % opt
par.opt_PStor_scale  = off ;
par.opt_alphaS       = off ; % opt
par.opt_gammaDNA	 = off ;
% O model parameters
par.opt_O2C_T = off ;
par.opt_rO2C  = off ;
% Si model parameters
par.opt_dsi   = off  ;
par.opt_at    = off ;
par.opt_bt    = off  ;
par.opt_aa    = off  ;
par.opt_bb    = off  ;
%
%-------------load data and set up parameters---------------------
SetUp ;                      

mmC = 12.0;         % molar mass of Carbon
mmP = 31;           % molar mass of P = 31 grams/mol

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

par.fname_diags = strcat(fname,'_diagnostics.mat') ;
% -------------------update initial guesses --------------
if isfile(par.fname)
    fprintf('loading initial guess on C and O from file: %s \n',par.fname)
    load(par.fname,'data')
end 

% -------------------update initial guesses --------------
if isfile(par.fnameload)
    fprintf('loading initial guess on C and O from file: %s \n',par.fnameload)
    load(par.fnameload,'data')
else
    fprintf('did not load a previous model solution from file par.fnameload \n')
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


%% 
%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPar(par)  ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

x0    = par.p0 ;
nip = length(x0);


% -------- initialize model state with optimal solution -----------
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

Xin.P = [par.DIP;par.POP;par.DOP;par.DOPl];
Xin.C = [par.DIC;par.POC;par.DOC;par.PIC;par.ALK;par.DOCl;par.DOCr]; %; par.pco2atm];
Xin.O2 = [par.O2];



%% -------------------set up timestepper -------------------------
spd = 60*60*24; 
spa  = 365*spd ;
mmC = 12.0;
mmP = 31.0;

M3dsurf = M3d ;             % make surface mask (ocn grid cells in contact with atm) 
M3dsurf(:,:,2:end) = 0 ;
Msurf = M3dsurf(iwet);
isrf = find(M3dsurf(iwet)) ;

M3dEZ = M3d ;             % make surface mask (ocn grid cells in contact with atm) 
M3dEZ(:,:,3:end) = 0 ;
iEZ = find(M3dEZ(iwet)) ;

M3dDeep = M3d;
M3dDeep(:,:,1:2) = 0 ;
iDeep = find(M3dDeep(iwet));

%%
% Default schedule
%                      1 hr    1 day   1 month  1 year   4 years   
dt_size =        spa./[365*24  365     12       1        0.25   ];  % step sizes in seconds
nsteps =              [ 24     364     50       20       30     ]; % number of steps for each size

%                      6 hr      4 day   6 hr      4 day.  1 month  1 year   4 years  
dt_size =        spa./[365*24/6  365/4   365*24/6  365/4     12       1        0.25   ];  % step sizes in seconds
nsteps =              [ 4        91        4       364/4     48       45       25     ]; % number of steps for each size

%                      1 hr      1 day   1 hr     1 day.  1 month  1 year   4 years  
dt_size =        spa./[365*24    365     365*24   365     12        1        0.25   ];  % step sizes in seconds
nsteps =              [ 240      355     240      355     120       90       25     ]; % number of steps for each size

%% test run
%dt_size =       spa./[365*24    365]; %    12    ]; %1    0.25];  % step sizes in seconds
%nsteps =              [ 240     355   ];  %24    ]; %50   25]; % number of steps for each size

Nstep_save = 10; % Number of steps between saving output

par.dt_schedule = dt_size ;
par.nsteps_schedule = nsteps ;

total_nsteps = sum(nsteps);
total_time_yr = sum(dt_size .* nsteps) / spa;
fprintf('Time-step Schedule: %d total steps over %.2f years\n', total_nsteps, total_time_yr);


% allocate output  
% diagnostics to save at every timestep: 
% pco2atm, volume integrated DIC, integrated PNPP, CNPP
tmp = zeros(1,total_nsteps);
Tout = zeros(1,total_nsteps);
pco2atmout = zeros(1,total_nsteps);
Diags.Tout = Tout ;
Diags.PNPP = tmp;
Diags.CNPP_nolabile = tmp;
Diags.pco2atm = pco2atmout ;
Diags.totalDIC = tmp; 
Diags.totalDICdeep = tmp;
Diags.totalDICsurf = tmp;
Diags.totalDOC = tmp;
Diags.totalDOCl = tmp;
Diags.totalDOCr = tmp;
Diags.totalPOC = tmp;
Diags.totalPIC = tmp;
Diags.totalDIP = tmp; 
Diags.totalDIPdeep = tmp;  
Diags.totalDIPsurf = tmp;
Diags.totalPOP = tmp; 
Diags.totalDOP = tmp; 
Diags.totalDOPl = tmp; 
Diags.totalPexport_darkremin = tmp;
Diags.totalCexport_darkremin = tmp;


% initialize output structure
if par.saveall
    OUT.P = zeros(length(Xin.P),total_nsteps); 
    OUT.C = zeros(length(Xin.C),total_nsteps); 
end



t0 = 0; 
  t = t0;
  Tout(1) = t;
  fprintf('Before time stepping method (Steady-state solutions in pre-industrial era)...\n');
  fprintf('...Time: %4.2f, AtmCO2: %4.2f uatm,  avgDIC: %7.6g mmol/m3\n', ...
    t,par.pco2atm,sum(par.DIC.*par.dVt(iwet))/sum(par.dVt(iwet))); % mean(par.DIC)

    totalDIC = sum(par.DIC.*par.dVt(iwet).*1e-3); % units = mol C
    totalCO2atm = par.pco2atm * par.Natm * 1e-6; % mol C
    fprintf('...Atm: Integrated total CO2  = %10.3e Pg C \n',totalCO2atm*12*1e-15);
    %fprintf('...Ocn: avgDIC = %7.4f mmol/m3 \n',sum(par.DIC.*par.dVt(iwet))/sum(par.dVt(iwet)));
    fprintf('...Ocn: Integrated total DIC  = %10.3e Pg C \n',totalDIC*12*1e-15); 
    fprintf('...Ocn: Integrated total POC  = %10.3e Pg C \n',sum(par.POC.*par.dVt(iwet).*1e-3)*12*1e-15);    
    fprintf('...Ocn: Integrated total DOC  = %10.3e Pg C \n',sum(par.DOC.*par.dVt(iwet).*1e-3)*12*1e-15);
    fprintf('...Ocn: Integrated total DOCl = %10.3e Pg C \n',sum(par.DOCl.*par.dVt(iwet).*1e-3)*12*1e-15);
    fprintf('...Ocn: Integrated total DOCr = %10.3e Pg C \n',sum(par.DOCr.*par.dVt(iwet).*1e-3)*12*1e-15);
    fprintf('...Ocn: Integrated total PIC  = %10.3e Pg C \n',sum(par.PIC.*par.dVt(iwet).*1e-3)*12*1e-15);

%% change npp to simulate fertilization.
% calculate NPP at steady state
    

    LAM = 0*M3d;
    for ji = 1 : par.nl
        LAM(:,:,ji) = (par.npp(:,:,ji).^par.beta).*par.Lambda(:,:,ji) ;
    end 

    L      = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
    par.L  = L;
    
    % organic P production
    G = M3d*0;
    G(iwet) = par.alpha*L*par.DIP;

    C2P3D = data.C2P;

    % organic C production
    % not including labile
    Cprod_pocdoc = G.*C2P3D;
    % par.Cnpp; % unit: (mmolC/m^3/s) 
    Cprod_docl = par.Cnpp - Cprod_pocdoc;

    tmp = sum(par.Cnpp(iwet).*dVt(iwet),'all','omitnan')*mmC*spd*365*1e-18;
    fprintf('Global satellite C NPP: %3.2f Pg C /yr \n\n',tmp);

    globalCprod_pocdoc = sum(Cprod_pocdoc(iwet).*dVt(iwet),'all','omitnan')*mmC*spd*365*1e-18;
    fprintf('Global POC and DOC production (G*C2P) excluding labile DOC: %3.2f Pg C /yr \n',globalCprod_pocdoc);

    tmp = sum(Cprod_docl(iwet).*dVt(iwet),'all','omitnan')*mmC*spd*365*1e-18;
    fprintf('Global labile DOC production (satellite CNPP - G*C2P): %3.2f Pg C /yr \n\n',tmp);
%%
    % % increase L, DIP uptake rate, by 10% everywhere
    % fprintf('10 percent increase L ... \n')
    % G = M3d*0;
    % G(iwet) = 1.1*par.alpha*L*par.DIP;
    % Cprod_pocdoc = G.*C2P3D;
    % globalCprod_pocdoc = sum(Cprod_pocdoc(iwet).*dVt(iwet),'all','omitnan')*mmC*spd*365*1e-18;
    % fprintf('Global POC and DOC production (G*C2P) excluding labile DOC: %3.2f Pg C /yr \n',globalCprod_pocdoc);

    % 2.6 Pg C increase. but biggest changes in subtropical gyres where L is
    % high, but this is not where fertilization would happen
    % Instead, make a mask to only increase L in southern ocean?
    % make an HNLC mask on 91x180x24 grid, only need surface 2 layers for production
    % make a Southern ocean mask
    MSKS.SO = M3d; % initialize mask of all ocean grid points   
    MSKS.SO(find(grd.yt>-45),:,:) = 0; % set everything north of -45 lat to zero


    % Equatorial Pacific HNLC
    % lat bounds -10:10
    % lon bounds 140:270 
    MSKS.EqP = M3d;

    MSKS.EqP(find(grd.yt>5),:,:) = 0;
    MSKS.EqP(find(grd.yt<-5),:,:) = 0;
    MSKS.EqP(:,find(grd.xt>270),:) = 0;
    MSKS.EqP(:,find(grd.xt<140),:) = 0;

    % N. Pacific high lat
    %lon: 140:220
    %lat bounds 45:65 N
    MSKS.NP = M3d;
    MSKS.NP(find(grd.yt>65),:,:) = 0;
    MSKS.NP(find(grd.yt<45),:,:) = 0;
    MSKS.NP(:,find(grd.xt<140),:) = 0;
    MSKS.NP(:,find(grd.xt>220),:) = 0;

    %figure; imagesc(MSKS.EqP(:,:,1)); axis xy; cb = colorbar; title('Equatorial Pacific HNLC mask')
    %figure; imagesc(MSKS.NP(:,:,1)); axis xy; cb = colorbar; title('North Pacific HNLC mask')
    %figure; imagesc(M3d(:,:,1)); axis xy; cb = colorbar; title('ocean mask')

    MSKS.HNLC = (MSKS.SO | MSKS.EqP | MSKS.NP);
    figure; 
    imAlpha = M3d(:,:,1); imagesc(MSKS.HNLC(:,:,1),'AlphaData',imAlpha); axis xy; cb = colorbar; title('All HNLC mask')

%% now increase L in only HNLC regions
    % increase L, DIP uptake rate, by 10% everywhere
    fprintf('50 percent increase L in HNLC ... \n')
    G = M3d*0;
    G(iwet) = par.alpha*L*par.DIP;
    % increase uptake rate in HNLC
    G(MSKS.HNLC) = 1.5*G(MSKS.HNLC);
    Cprod_pocdoc = G.*C2P3D;
    globalCprod_pocdoc = sum(Cprod_pocdoc(iwet).*dVt(iwet),'all','omitnan')*mmC*spd*365*1e-18;
    fprintf('Global POC and DOC production (G*C2P) excluding labile DOC: %3.2f Pg C /yr \n',globalCprod_pocdoc);

    % 27.92 Pg C/yr (from 25.96 Pg C/yr)

    % Implement as change to par.Lambda, since that is created in SetUp and used in Peqn
    % or implemenet as change to par.npp, but thats a bit harder since par.npp is raised to pawer of beta
    % after running setup
    % par.Lambda(jj,ii,1) = 1./(1e-6+DIP_obs(jj,ii,1)) ;      % unit: [1/(mmolP/m^3)]
    % par.Lambda(jj,ii,2) = 1./(1e-6+DIP_obs(jj,ii,2)) ;
    % in Peqn, LAM(:,:,ji) = (npp(:,:,ji).^beta).*Lambda(:,:,ji) ; L = d0(LAM(iwet))

    % for first year of timesteps, increase Lambda in HNLC regions by 50%
    % par.Lambda(MSKS.HNLC) = 1.5*par.Lambda(MSKS.HNLC);

    % after 1 year of timesteps, reset Lambda to original value
    % par.Lambda(MSKS.HNLC) = (1/1.5)*par.Lambda(MSKS.HNLC);


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
fprintf('Solve eqPcycle...\n')
[par, P] = eqPcycle(x0, par) ;
    par.DIP = P(1:nwet);
    par.POP = P(1*nwet+1:2*nwet);
    par.DOP = P(2*nwet+1:3*nwet);
    par.DOPl= P(3*nwet+1:4*nwet);

    % check if P(1:4*nwet) equals Xin.P(1:4*nwet)
    fprintf('Comparing P to Xin.P...\n');
    if numel(P) >= 4*nwet && numel(Xin.P) >= 4*nwet
        A = P(1:4*nwet);
        B = Xin.P(1:4*nwet);
        if isequal(A,B)
        fprintf('P(1:4*nwet) exactly equals Xin.P(1:4*nwet)\n');
        else
        D = A - B;
        maxabs = max(abs(D));
        maxrel = max(abs(D) ./ (abs(B) + eps));
        fprintf('P vs Xin.P: max abs diff = %.3e, max rel diff = %.3e\n', maxabs, maxrel);
        end
    else
        warning('P or Xin.P does not contain 4*nwet elements (have %d and %d)', numel(P), numel(Xin.P));
    end

    % check that eqCcycleAtm solution matches eqCcycle_v2 solution for DIC, POC, DOC, PIC, ALK, DOCl, DOCr
    try
        % call eqCcycle_v2 (signature used in this repo: [par,C,...] = eqCcycle_v2(x,par))
        fprintf('Solve eqCcycle_v2...\n');
        if isfile('../output/PNAS2025_transient/transient_test_steadystate_Conly_1h_Atm_from_reoptNature_with_dop_GM15_npp1_CTL_He_PC_parXXXX.mat')
            fprintf('Loading presaved par file for eqCcycle_v2...\n');
            %Code_withCellModel/output/PNAS2025_transient/transient_test_steadystate_Conly_1h_noAtm_from_reoptNature_with_dop_GM15_npp1_CTL_He_PC_par.mat
            %../output/PNAS2025_transient/transient_test_steadystate_Conly_1h_Atm_from_reoptNature_with_dop_GM15_npp1_CTL_He_PC_par.mat ...
            par_tmp = load('../output/PNAS2025_transient/transient_test_steadystate_Conly_1h_Atm_from_reoptNature_with_dop_GM15_npp1_CTL_He_PC_par.mat','par');
            par.DIC = par_tmp.par.DIC;
            par.POC = par_tmp.par.POC;
            par.DOC = par_tmp.par.DOC;
            par.PIC = par_tmp.par.PIC;
            par.ALK = par_tmp.par.ALK;
            par.DOCl= par_tmp.par.DOCl;
            par.DOCr= par_tmp.par.DOCr;
            GC  = [par.DIC; par.POC; par.DOC; par.PIC; ...
                   par.ALK; par.DOCl; par.DOCr];
        else
        GC  = [par.DIC; par.POC; par.DOC; par.PIC; ...
           par.ALK; par.DOCl; par.DOCr];
        
        [par_tmp, C_v2] = eqCcycle_v2(x0, par);
        
        fprintf('...eqCcycle_v2  done.\n');
        par_tmp.DIC = C_v2(1:nwet);
        par_tmp.POC = C_v2(1*nwet+1:2*nwet);
        par_tmp.DOC = C_v2(2*nwet+1:3*nwet);
        par_tmp.PIC = C_v2(3*nwet+1:4*nwet);
        par_tmp.ALK = C_v2(4*nwet+1:5*nwet);
        par_tmp.DOCl= C_v2(5*nwet+1:6*nwet);
        par_tmp.DOCr= C_v2(6*nwet+1:7*nwet);
        %par_tmp.pco2atm= C_v2(7*nwet+1);

        % check if C_v2(1:7*nwet) equals Xin.C(1:7*nwet)
        fprintf('Comparing C_v2 to Xin.C...\n');
        if numel(C_v2) >= 7*nwet && numel(Xin.C) >= 7*nwet
            A = C_v2(1:7*nwet);
            B = Xin.C(1:7*nwet);
            if isequal(A,B)
            fprintf('C_v2(1:7*nwet) exactly equals Xin.C(1:7*nwet)\n');
            else
            D = A - B;
            maxabs = max(abs(D));
            maxrel = max(abs(D) ./ (abs(B) + eps));
            fprintf('C_v2 vs Xin.C: max abs diff = %.3e, max rel diff = %.3e\n', maxabs, maxrel);
            end
        else
            warning('C_v2 or Xin.C does not contain 7*nwet elements (have %d and %d)', numel(C_v2), numel(Xin.C));
        end
        
        par = par_tmp; 
    end

    catch ME
        warning('Could not run eqCcycle_v2: %s', ME.message);
        C_v2 = [];
    end



% fprintf('Solve eqCcycleAtm...\n');
% GC  = [par.DIC; par.POC; par.DOC; par.PIC; ...
%            par.ALK; par.DOCl; par.DOCr; par.pco2atm];
% %[par, C, Cx, Cxx] = eqCcycle_v2(x, par)
% [par, C] = eqCcycleAtm(x0, par);
%     par.DIC    = C(1:nwet);
%     par.POC    = C(1*nwet+1:2*nwet);
%     par.DOC    = C(2*nwet+1:3*nwet);
%     par.PIC    = C(3*nwet+1:4*nwet);
%     par.ALK    = C(4*nwet+1:5*nwet);
%     par.DOCl   = C(5*nwet+1:6*nwet);
%     par.DOCr   = C(6*nwet+1:7*nwet);
%     par.pco2atm= C(7*nwet+1);

    fprintf('New eqCcycle Solution \n')
    fprintf('...Time: %4.2f, AtmCO2: %4.2f uatm,  avgDIC: %7.6g mmol/m3\n', ...
    t,par.pco2atm,sum(par.DIC.*par.dVt(iwet))/sum(par.dVt(iwet))); % mean(par.DIC)
    totalDIC = sum(par.DIC.*par.dVt(iwet).*1e-3); % units = mol C
    totalCO2atm = par.pco2atm * par.Natm * 1e-6; % mol C
    fprintf('...Atm: Integrated total CO2  = %10.3e Pg C \n',totalCO2atm*12*1e-15);
    fprintf('...Ocn: Integrated total DIC  = %10.3e Pg C \n',totalDIC*12*1e-15); 
    fprintf('...Ocn: Integrated total POC  = %10.3e Pg C \n',sum(par.POC.*par.dVt(iwet).*1e-3)*12*1e-15);    
    fprintf('...Ocn: Integrated total DOC  = %10.3e Pg C \n',sum(par.DOC.*par.dVt(iwet).*1e-3)*12*1e-15);
    fprintf('...Ocn: Integrated total DOCl = %10.3e Pg C \n',sum(par.DOCl.*par.dVt(iwet).*1e-3)*12*1e-15);
    fprintf('...Ocn: Integrated total DOCr = %10.3e Pg C \n',sum(par.DOCr.*par.dVt(iwet).*1e-3)*12*1e-15);
    fprintf('...Ocn: Integrated total PIC  = %10.3e Pg C \n',sum(par.PIC.*par.dVt(iwet).*1e-3)*12*1e-15);

    fprintf('\n');
    %save('temp_par_before_eqCcycleAtm.mat','par');
    fprintf('save temporary par, after eqPcycle and eqCcyle_v2 to file: %s ...\n',par.fxpar);
    save(par.fxpar,'par');

%     if ~isempty(C_v2)
%         species = {'DIC','POC','DOC','PIC','ALK','DOCl','DOCr'};
%         tol_rel = 1e-8;
%         tol_abs = 1e-12;
%         ok_all = true;
%         for k = 1:7
%             i1 = (k-1)*nwet + 1;
%             i2 = k*nwet;
%             A = C(i1:i2);
%             B = C_v2(i1:i2);
%             if numel(A) ~= numel(B)
%                 warning('Size mismatch for %s: %d vs %d', species{k}, numel(A), numel(B));
%                 ok_all = false;
%                 continue
%             end
%             D = A - B;
%             maxabs = max(abs(D));
%             maxrel = max(abs(D)./(abs(B) + eps));
%             fprintf('%s: max abs diff = %.3e , max rel diff = %.3e\n', species{k}, maxabs, maxrel);
%             if maxabs > tol_abs && maxrel > tol_rel
%                 warning('%s difference exceeds tolerances (abs>%.1e and rel>%.1e)', species{k}, tol_abs, tol_rel);
%                 ok_all = false;
%             end
%         end
%         % compare atmospheric pCO2 if present
%         if numel(C) >= 7*nwet+1 && numel(C_v2) >= 7*nwet+1
%             ap = C(7*nwet+1);
%             bp = C_v2(7*nwet+1);
%             pd = ap - bp;
%             fprintf('pco2atm: abs diff = %.3e , rel diff = %.3e\n', abs(pd), abs(pd)./(abs(bp)+eps));
%             if abs(pd) > tol_abs && abs(pd)./(abs(bp)+eps) > tol_rel
%                 warning('pco2atm difference exceeds tolerances');
%                 ok_all = false;
%             end
%         end
%         if ok_all
%             fprintf('eqCcycleAtm and eqCcycle_v2 match within tolerances.\n');
%         else
%             fprintf('Differences found between eqCcycleAtm and eqCcycle_v2.\n');
%             keyboard;
%         end
%     end
%%
% check the data and par for steady-state drive the time stepping method 

% change npp to simulate fertilization.

if par.equalbkCP_flag ==true
    %% reset parameters to make remin rates equal for P and C
    % par.fxhatload = '../output/test_equalbkCP_reoptNature_with_dop_GM15_npp1_CTL_He_xhat.mat
    fprintf('Setting P remin parameters (kdP,Q10P,bP,bP_T) equal to C remin parameters \n')
    par.kdP = par.kdC ; % set P remin rate equal to C remin rate
    par.Q10P = par.Q10C ; % set P Q10 equal to C Q10
    par.bP = par.bC ; % set P b equal to C b
    par.bP_T = par.bC_T ; % set P b_T equal to C b_T
end

%% set up parameters for diagnostic calculations
    nl = 2; % number of layers in euphotic zone 
    tf      = (par.vT - 30)/10 ;
    kP      =  par.kdP * par.Q10P.^tf ;
    kP3d    = M3d+nan;
    kP3d(iwet)      = kP ;
    kC      = d0(par.kdC * par.Q10C .^ tf) ;
    kC3d    = M3d+nan;
    kC3d(iwet)      = par.kdC * par.Q10C .^ tf ;
    kappa_r = par.kru*par.UM + par.krd*par.DM ;
    eta     = par.etau*par.WM ;
    par.kappa_r = kappa_r ; 
    par.eta     = eta ;
    par.kP = kP; 
    par.kC = kC; 
 


% set up iteration counter 
t = t0;
global_step = 0; % global step counter for saving
for dt_idx = 1:length(dt_size)
    dt = dt_size(dt_idx);
    n_substeps = nsteps(dt_idx);
    % load previous iteration solution
    Xin.P = [par.DIP;par.POP;par.DOP;par.DOPl];
    if numel(Xin.C) >= 7*nwet+1
        Xin.C = [par.DIC;par.POC;par.DOC;par.PIC;par.ALK;par.DOCl;par.DOCr;par.pco2atm];
    else
        Xin.C = [par.DIC;par.POC;par.DOC;par.PIC;par.ALK;par.DOCl;par.DOCr];
    end
    Xin.O2 = [par.O2];

    % Reset Lambda to simulate fertilization
    par.Lambda = M3d*0 ;
    [nx,ny,nz] = size(M3d) ;
    for jj = 1 : nx
        for ii = 1 : ny     
            par.Lambda(jj,ii,1) = 1./(1e-6+DIP_obs(jj,ii,1)) ;      % unit: [1/(mmolP/m^3)]
            par.Lambda(jj,ii,2) = 1./(1e-6+DIP_obs(jj,ii,2)) ;
        end
    end

    %increase lambda by 50% in HNLC regions, for first year (2 )
    if dt_idx < 3 & par.fertilize_flag == true;                  
        % if fertilization == on, increase P uptake by enough to draw down surface DIP in S.O.
        par.Lambda(MSKS.HNLC) = 1.5*par.Lambda(MSKS.HNLC);
        fprintf('...Increase Lambda in HNLCs by 50 percent from steady state model\n')
    else
        fprintf('...Same Lambda as steady state model\n')
    end

    % set up time stepper for P and C (only need to factor trapezoid matrix once for each step size dt)
    % run Peqn
    fprintf('\ndt = %.1f \n', dt)
    fprintf('...run Peqn \n')
    [f_P,J_P,par] = Peqn(Xin.P, par); 
    % Evaluate RHS and Jacobian at current state (X, t)
    % Build trapezoidal matrices
    I_P  = speye(numel(Xin.P(:)));
    A_P  = I_P + 0.5*dt*J_P;
    B_P  = I_P - 0.5*dt*J_P;
    fprintf('...factor the big matrix... \n')
    tic
    A_Pfactored = mfactor(A_P); 
    toc

    % run CeqnAtm
    % fprintf('...run CeqnAtm \n')
    % [f_C,J_C,par] = CeqnAtm(Xin.C, par);
    fprintf('...run Ceqn_v2 \n')
    [f_C,J_C,par] = Ceqn_v2(Xin.C, par);
    % Evaluate RHS and Jacobian at current state (X, t)
    % Build trapezoidal matrices
    I_C  = speye(numel(Xin.C(:)));
    A_C  = I_C + 0.5*dt*J_C;
    B_C  = I_C - 0.5*dt*J_C;
    fprintf('...factor the big matrix... \n')
    tic
    A_Cfactored = mfactor(A_C);
    toc

    % start iterating through n_substeps for each dt
    tic
    for ii = 1:n_substeps 
        global_step = global_step + 1; 
        %fprintf('%i. ',ii); 
        t = t+dt;
        Tout(global_step) = t;

        % load previous iteration solution
        Xin.P = [par.DIP;par.POP;par.DOP;par.DOPl];
        if numel(Xin.C) >= 7*nwet+1
            Xin.C = [par.DIC;par.POC;par.DOC;par.PIC;par.ALK;par.DOCl;par.DOCr;par.pco2atm];
        else
            Xin.C = [par.DIC;par.POC;par.DOC;par.PIC;par.ALK;par.DOCl;par.DOCr];
        end
        Xin.O2 = [par.O2];

        % run Peqn
        [f_P,J_P,par] = Peqn(Xin.P, par); 
        % Right-hand side
        rhs_P = B_P*Xin.P - dt*f_P;
        Xout.P  = mfactor(A_Pfactored, rhs_P);        

        % update par with Peqn solution
        par.DIP = Xout.P(1:nwet);
        par.POP = Xout.P(1*nwet+1:2*nwet);
        par.DOP = Xout.P(2*nwet+1:3*nwet);
        par.DOPl= Xout.P(3*nwet+1:4*nwet);

        % % run CeqnAtm
        %[f_C,J_C,par] = CeqnAtm(Xin.C, par);
        [f_C,J_C,par] = Ceqn_v2(Xin.C, par);
        % Right-hand side
        %C=mfactor(FAC,BC*C-dt*fC);
        rhs_C = B_C*Xin.C - dt*f_C;
        Xout.C  = mfactor(A_Cfactored, rhs_C);
        % update par with CeqnAtm solution
        par.DIC    = Xout.C(1:nwet);
        par.POC    = Xout.C(1*nwet+1:2*nwet);
        par.DOC    = Xout.C(2*nwet+1:3*nwet);
        par.PIC    = Xout.C(3*nwet+1:4*nwet);
        par.ALK    = Xout.C(4*nwet+1:5*nwet);
        par.DOCl   = Xout.C(5*nwet+1:6*nwet);
        par.DOCr   = Xout.C(6*nwet+1:7*nwet);
        if numel(Xout.C) >= 7*nwet+1
            par.pco2atm= Xout.C(7*nwet+1);
        end

        

        % save output at each time step
        if par.saveall
            OUT.P(:,global_step) = Xout.P;
            OUT.C(:,global_step) = Xout.C;
        else

        end

        % if ~par.saveall
        %     if mod(global_step,Nstep_save) == 0
        %     save_step = save_step +1;
        %     OUT.P{save_step} = Xout.P;
        %     outname = sprintf('%s%sTstep_%i_%.3fyr.mat',output_dir,VerName,global_step,t/spa);
        %     fprintf('saving model step to file: %s \n',outname)
        %     save(outname, 'Xout','Tout','t','global_step','-v7.3')
        %     end
        % end
        
        Diags.Tout(global_step) = Tout(global_step);
        Diags.pco2atm(global_step) = par.pco2atm;        

        % calculate NPP
        PNPP = par.alpha*par.L*par.DIP;
        CNPP_nolabile = PNPP.*par.C2P;

        Diags.PNPP(global_step) = sum(PNPP.*dVt(iwet),'all','omitnan')*mmP*spd*365*1e-18 ; 
        Diags.CNPP_nolabile(global_step) = sum(CNPP_nolabile.*dVt(iwet),'all','omitnan')*mmC*spd*365*1e-18 ; 

        % calculate integrated DIC inventory in the ocean
        totalDIC = sum(par.DIC.*par.dVt(iwet).*1e-3)*12*1e-15; % units = Pg C

        Diags.totalDIC(global_step) = totalDIC;
        Diags.totalDICsurf(global_step) = sum(par.DIC(iEZ).*dVt(iwet(iEZ)))*mmC*1e-3*1e-15; %Pg DIC
        Diags.totalDICdeep(global_step) = sum(par.DIC(iDeep).*dVt(iwet(iDeep)))*mmC*1e-3*1e-15; %Pg DIC
     
        % calculate total POP in euphotic zone
        %totalPOPsurf = sum(par.POP(iEZ).*dVt(iwet(iEZ)))*mmP*1e-3*1e-15; %Pg POP

        % calculate total POC in euphotic zone
        %Diags.totalPOCsurf(global_step) = sum(par.POC(iEZ).*dVt(iwet(iEZ)))*mmC*1e-3*1e-15; %Pg POC

        Diags.totalDOC(global_step) = sum(par.DOC.*dVt(iwet))*mmC*1e-3*1e-15; %Pg DOC
        Diags.totalDOCl(global_step) = sum(par.DOCl.*dVt(iwet))*mmC*1e-3*1e-15; %Pg DOC
        Diags.totalDOCr(global_step) = sum(par.DOCr.*dVt(iwet))*mmC*1e-3*1e-15; %Pg DOC
        Diags.totalPOC(global_step) = sum(par.POC.*dVt(iwet))*mmC*1e-3*1e-15; %Pg POC
        Diags.totalPIC(global_step) = sum(par.PIC.*dVt(iwet))*mmC*1e-3*1e-15; %Pg PIC

        % calculate surface DIP
        Diags.totalDIP(global_step) = sum(par.DIP.*dVt(iwet))*mmP*1e-3*1e-15 ;
        Diags.totalDIPsurf(global_step) = sum(par.DIP(iEZ).*dVt(iwet(iEZ)))*mmP*1e-3*1e-15; %Pg DIP
        Diags.totalDIPdeep(global_step) = sum(par.DIP(iDeep).*dVt(iwet(iDeep)))*mmP*1e-3*1e-15; %Pg DIP
        Diags.totalPOP(global_step) = sum(par.POP.*dVt(iwet))*mmP*1e-3*1e-15; %Pg POP
        Diags.totalDOP(global_step) = sum(par.DOP.*dVt(iwet))*mmP*1e-3*1e-15; %Pg DOP
        Diags.totalDOPl(global_step) = sum(par.DOPl.*dVt(iwet))*mmP*1e-3*1e-15; %Pg DOPl
    
        % store 3D fields
        model = struct();
        model.DIP = M3d+nan ; model.POP = M3d+nan ; model.DOP = M3d+nan ; model.DOPl = M3d+nan ;
        model.DIP(iwet) = par.DIP ;
        model.POP(iwet) = par.POP ;
        model.DOP(iwet) = par.DOP ;  
        model.DOPl(iwet) = par.DOPl ; 

        model.DIC = M3d+nan ; model.POC = M3d+nan ; model.DOC = M3d+nan ; 
        model.PIC = M3d+nan ; model.ALK = M3d+nan ; model.DOCr = M3d+nan ; model.DOCl = M3d+nan ;
        model.DIC(iwet) = par.DIC ;
        model.POC(iwet) = par.POC ;
        model.DOC(iwet) = par.DOC ;
        model.PIC(iwet) = par.PIC ;
        model.ALK(iwet) = par.ALK ;
        model.DOCr(iwet) = par.DOCr ;
        model.DOCl(iwet) = par.DOCl ;

        % calculate P export (water-column integrated remineralization)
        % Total organic P remineralization
        orgPremin = (kP3d.*model.DOP + par.kappa_p*model.POP + par.kappa_l*model.DOPl);         % [mmol P/m^3/sec]
        Pexp = orgPremin(:,:,nl+1:end).*dVt(:,:,nl+1:end);                              % [mmol P/s]
        Pexpint = nansum(Pexp,3);                                                       % vertical sum [mmol P/s]
        Diags.totalPexport_darkremin(global_step) = nansum(Pexpint(:))*mmP*spa*1e-18;                                  % global sum [Pg P/yr]
        %fprintf('Model globally-integrated TOP  remineralization below the Euphotic zone is %3.4f Pg P /yr \n\n',Diags.totalPexport_darkremin);
        %Diags.TOPexp_wcremin{global_step} = sum(orgPremin(:,:,nl+1:end).*grd.DZT3d(:,:,3:end),3,'omitnan')*mmP*spd; %[mg P/m^2/day] vertical sum
        

        % Total organic C remineralization
        orgCremin = M3d+nan;
        orgCremin(iwet) = par.kappa_p*model.POC(iwet) + eta*(kC*model.DOC(iwet)) + par.kappa_l*model.DOCl(iwet) +par.kappa_r*model.DOCr(iwet); % [mmol C/m^3/s]
        Cexp = orgCremin(:,:,nl+1:end).*dVt(:,:,nl+1:end);                  % [mmol C/sec]
        Cexpint = sum(Cexp,3,'omitnan');                                    % vertical sum [mmol C/s]
        Diags.totalCexport_darkremin(global_step) = sum(Cexpint(:),'all','omitnan')*mmC*spa*1e-18;     % global sum [Pg C/yr]
        %fprintf('Model globally-integrated TOC  remineralization below the Euphotic zone is %3.3f Pg C /yr \n\n',Diags.totalCexport_darkremin);


        % concise diagnostic printout
        fprintf('t=%6.1f d | pCO2=%6.2f ppm | totalDIC=%8.3f PgC | DIP_surf=%6.3e PgP | PNPP=%6.3f PgP/yr | CNPP_nolab=%6.3f PgC/yr |  orgCexp=%6.3e Pg C/yr\n', ...
            t/spd, par.pco2atm, totalDIC, Diags.totalDIPsurf(global_step), Diags.PNPP(global_step), Diags.CNPP_nolabile(global_step),  Diags.totalCexport_darkremin(global_step));

        % fprintf('t= %10.2f, Ocn: totalDIC = %11.4e Pg C ; Atm pCO2 = %10.3f ppm \n',t/spd, totalDIC, par.pco2atm)

    end
    toc
    %sprintf('%sdt_%.0fd.mat', par.VerName,dt/spd) 
end

%[par.VerName 'diagnostics.mat']
fprintf('saving diagnostics to file: %s \n',par.fname_diags)
save(par.fname_diags,"Diags","Tout","dt_size","nsteps")

if par.saveall
    fprintf('saving model solution to file: %s \n',par.fname)
    save(par.fname, 'Tout','OUT','-v7.3');
end



%  if exist(par.fname, 'file')
%     reply = input(sprintf('WARNING: File ( %s ) already exists. \nDo you want to overwrite this file? Y/N: ', par.fname), 's');
%     if strcmpi(reply, 'Y')
%         fprintf('Overwriting File... \n');
%         fprintf('saving model solution to file: %s \n',par.fname)
%         save(par.fname, 'data')
%     else
%         fprintf('Execution stopped by User.\n');
%         fprintf('--------------------------\n\n');
%         return;
%     end
% else
%     fprintf('saving model solution to file: %s \n',par.fname)
%     save(par.fname, 'data')
% end

%save output after 1 year of fertilization. 
% start new timestep to continue without fertilization



% save model solution to data file
        % DIC  = M3d+nan ;  DIC(iwet)  = par.DIC ;
        % POC  = M3d+nan ;  POC(iwet)  = par.POC ;
        % DOC  = M3d+nan ;  DOC(iwet)  = par.DOC ;
        % PIC  = M3d+nan ;  PIC(iwet)  = par.PIC ;
        % ALK  = M3d+nan ;  ALK(iwet)  = par.ALK ;
        % DOCl = M3d+nan ;  DOCl(iwet) = par.DOCl ;
        % DOCr = M3d+nan ;  DOCr(iwet) = par.DOCr ;
        % 
        % data.DIC  = DIC  ;  data.POC  = POC ;
        % data.DOC  = DOC  ;  data.PIC  = PIC ;
        % data.ALK  = ALK  ;  data.DOCr = DOCr ;
        % data.DOCl = DOCl ;  data.pco2atm = par.pco2atm ;
        % fprintf('saving model solution to file: %s \n',par.fname)
        % save(par.fname, 'data')


