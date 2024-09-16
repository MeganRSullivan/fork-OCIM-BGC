% driver_optC_Cell.m
%   Driver script to run biogeochemical model, without optimizing parameters.
%   This run only solves the phosphorus cycle, using predetermined parameters
% ------------------------------------------------------------------------
%clc; clear all; close all
global iter
iter = 0 ;
on   = true  ; off  = false ;
format short

% --- addpath to model code -----
%addpath('../src/')
addpath('../src_reoptNature/')

% test1_eqPcycle_with_DOPl_gamma1pct_from_reoptNature_with_dop_GM15_npp1

VerName = 'test1_eqPcycle_with_DOPl_gamma10pct_from_reoptNature_with_dop_GM15_npp1_'; 		% optional version name. leave as an empty character array
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
par.Cmodel  = off ; 
par.Omodel  = off ; 
par.Simodel = off ;
par.Cisotope  = off  ;
par.dynamicP = off ; % if on, cell model uses modeled DIP. if off, cell model uses WOA observed DIP field.
par.LoadOpt = on ; % if load optimial parameters. 
% to load parameter values from a run with a different name.
%par.fxhatload = '../../output/optPonly_CTL_He_P_xhat.mat';
%par.fxhatload = '/DFS-L/DATA/primeau/hojons1/Nature2023_BGC_reoptimized/src_Nature_parameter_Megan/MSK91/CTL_He_PCO_Gamma0_kl12h_O5_POC2DIC_GM15_Nowicki_npp1_aveTeu_diffSig_O2C_uniEta_DICrmAnthro_2L_Pnormal_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00_xhat.mat' 
par.fxhatload = '../output/reoptNature_with_dop_GM15_npp1_CTL_He_PCO_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00_xhat.mat'

% to use different model output for initial CO guess. 
%par.fnameload = '/DFS-L/DATA/primeau/hojons1/Nature2023_BGC_reoptimized/src_Nature_parameter_Megan/MSK91/CTL_He_PCO_Gamma0_kl12h_O5_POC2DIC_GM15_Nowicki_npp1_aveTeu_diffSig_O2C_uniEta_DICrmAnthro_2L_Pnormal_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00.mat' ;
par.fnameload = '../output/reoptNature_with_dop_GM15_npp1_CTL_He_PCO_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00.mat'

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
output_dir = sprintf('../output/'); 

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
if isfile(par.fnameload)
    fprintf('loading initial guess on C and O from file: %s \n',par.fnameload)
    load(par.fnameload)
end 

% -------------------update initial guesses --------------
if isfile(par.fname)
    fprintf('loading initial guess on C and O from file: %s \n',par.fname)
    load(par.fname)
end 

%---------------- inital guesses on C and O ---------------
if par.Cmodel == on 
    GC  = [data.DIC(iwet); data.POC(iwet); data.DOC(iwet); data.PIC(iwet); ...
           data.ALK(iwet); data.DOC(iwet); data.DOC(iwet);];
    GC  = real(GC) + 1e-6*randn(7*nwet,1) ;
end 
if par.Omodel == on 
    GO  = real(data.O2(iwet)) + 1e-9*randn(par.nwet,1);
end 

%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPar(par)  ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

%-------------------set up fminunc -------------------------
x0    = par.p0 ;
myfun = @(x) neglogpost(x, par);
objfuntolerance = 5e-11; %5e-12;
options = optimoptions(@fminunc                  , ...
                       'Algorithm','trust-region', ...
                       'GradObj','on'            , ...
                       'Hessian','on'            , ...
                       'Display','iter'          , ...
                       'MaxFunEvals',2000        , ...
                       'MaxIter',2000            , ...
                       'TolX', objfuntolerance   , ...     % 기존은 -7. decreasing
                       'TolFun',objfuntolerance  , ...     % 기존은 -7. decreasing
                       'DerivativeCheck','off'   , ...
                       'FinDiffType','central'   , ...
                       'PrecondBandWidth',Inf)   ;
%
nip = length(x0);
if (par.optim)
    % save SetUp fields
    fprintf('saving initial SetUp par structure to file: %s \n',par.fxpar)
    if exist(par.fxpar, 'file')
        reply = input(sprintf('WARNING: File ( %s ) already exists. \nDo you want to overwrite this file? Y/N: ', par.fxpar), 's');
        if strcmpi(reply, 'Y')
            fprintf('Overwriting File... \n');
            save(par.fxpar, 'par', '-v7.3');
        else
            fprintf('Execution stopped by User.\n');
            fprintf('--------------------------\n\n');
            return;
        end
    else
        save(par.fxpar, 'par', '-v7.3');
    end
    % optimize parameters
    [xsol,fval,exitflag] = fminunc(myfun,x0,options);
    fprintf('objective function tolerance = %5.1e \n',objfuntolerance);
    fprintf('----fminunc complete----\n')
    [f,fx,fxx,data,xhat] = neglogpost(xsol,par);
    fprintf('----neglogpost solved for final parameter values----\n')
    xhat.pindx = par.pindx;
    xhat.f   = f   ;
    xhat.fx  = fx  ;
    xhat.fxx = fxx ;
    % save results 
    fprintf('saving optimized parameters to file: %s \n',fxhat)
    fprintf('saving model solution to file: %s \n',par.fname)
    save(fxhat, 'xhat')
    save(par.fname, 'data')
else
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
    par.DIP  = DIP(iwet) ;
    data.DIP = DIP ; data.POP  = POP  ;
    data.DOP = DOP ; data.DOPl = DOPl ;
    % DIP & DOP error
    DOP = DOP + DOPl; % sum of semilabile and labile DOP ;
    eip = DIP(iwet(idip)) - par.po4raw(iwet(idip)) ;
    eop = DOP(iwet(idop)) - par.dopraw(iwet(idop)) ;
    f  = f + 0.5*(eip.'*Wip*eip) + 0.5*(eop.'*Wop*eop); 
    f_components.DIP = 0.5*(eip.'*Wip*eip);
    f_components.DOP = 0.5*(eop.'*Wop*eop); 

    data.f = f;
    data.f_components = f_components;

    fprintf('current objective function value is: %3.3e \n\n',f) 
    fprintf('current objective function value for fit to DIP is %3.3e \n',f_components.DIP) 
    fprintf('current objective function value for fit to DOP is %3.3e \n',f_components.DOP) 
    
    %%%%%%%%%%%%%%%%%%   End Solve P    %%%%%%%%%%%%%%%%%%%%

    %% note: skipping save for testing
    % if exist(par.fname, 'file')
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
        fprintf('saving model solution to file: %s \n',par.fname)
        save(par.fname, 'data')
    % end
end
fprintf('-------------- end! ---------------\n');
