function [RHS,J,par] = Peqn(X, par)    
%% unpack some useful stuff
    on = true; off = false;
    grd   = par.grd   ;
    M3d   = par.M3d   ;
    TRdiv = par.TRdiv ;
    iwet  = par.iwet  ;
    nwet  = par.nwet  ;
    dVt   = par.dVt   ;
    I     = par.I     ;
    Natm    = 1.773e20  ; %  molar volume of atmosphere
	vw    = dVt(iwet)./Natm ;

    M3dsurf = M3d ;             % make surface mask (ocn grid cells in contact with atm) 
	M3dsurf(:,:,2:end) = 0 ;
	Msurf = M3dsurf(iwet);
	isrf = find(M3dsurf(iwet)) ;
    
    % unpack previous iteration solution
    DIP  = X(0*nwet+1:1*nwet) ; 
    POP  = X(1*nwet+1:2*nwet) ; 
    DOP  = X(2*nwet+1:3*nwet) ; 
    DOPl = X(3*nwet+1:4*nwet) ; 

    % unpack the parameters
    sigP  = par.sigP;
    Q10P = par.Q10P;
    kdP  = par.kdP;
    bP_T = par.bP_T;
    bP  = par.bP;
    alpha  = par.alpha;
    beta  = par.beta;

    % fixed parameters
    vT      = par.vT ;
    DIPbar  = M3d(iwet)*par.DIPbar ;  % gobal arerage PO4 conc.[mmol m^-3];
    gamma   = par.gamma   ; % fraction goes to lDOM ;
    kappa_g = par.kappa_g ; % PO4 geological restore const.[s^-1] ;
    kappa_p = par.kappa_p ; % POP solubilization rate constant ;
    kappa_l = par.kappa_l ; % labile DOM remineralization rate [s^-1];

    npp     = par.npp ;     % net primary production
    tf      = (vT - 30)/10 ;
    kP      =  kdP * Q10P.^tf ;
    % build part of the biological DIP uptake operator
    Lambda = par.Lambda;
    LAM    = 0*M3d;

    % change NaN value to zero in npp and Lambda
    % ---> it happens due to DIP_obs has nan values
    % in Weilei's code, smoothit functon remove the nan values of PO4.
    Lambda(isnan(Lambda(:))) = 0;
    npp(isnan(npp(:))) = 0;

    for ji = 1 : par.nl
        LAM(:,:,ji) = (npp(:,:,ji).^beta).*Lambda(:,:,ji) ;
    end 

    L      = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
    par.L  = L;
   
    % particle flux
    PFD  = buildPFD(par,'POP');

    % Set up system of equations
    %FDIP
    eq1 = (TRdiv + alpha*L + kappa_g*I)*DIP - (kappa_p*I)*POP - d0(kP)*DOP -(kappa_l*I)*DOPl - kappa_g*DIPbar ;

    %FPOP
    eq2 = (-(1-sigP-gamma)*alpha*L)*DIP + (PFD + kappa_p*I)*POP ;

    %FDOP
    eq3 = (-sigP*alpha*L)*DIP + (TRdiv + d0(kP))*DOP ;

    %FDOPL
    eq4 = (-gamma*alpha*L)*DIP + (TRdiv + kappa_l*I)*DOPl ; 

    F = [eq1; eq2; eq3; eq4];


    % build Jacobian equations.
    % column 1 dF/dDIP
    Jp{1,1} = TRdiv + alpha*L + kappa_g*I ;
    Jp{2,1} = -(1-sigP-gamma)*alpha*L ;
    Jp{3,1} = -sigP*alpha*L  ;
    Jp{4,1} = -gamma*alpha*L  ;

    % column 2 dF/dPOP
    Jp{1,2} = -kappa_p*I      ;
    Jp{2,2} = PFD + kappa_p*I ;
    Jp{3,2} = 0*I             ;
    Jp{4,2} = 0*I             ;

    % column 3 dF/dDOP
    Jp{1,3} = -d0(kP)         ;
    Jp{2,3} = 0*I ;
    Jp{3,3} = TRdiv + d0(kP)  ;
    Jp{4,3} = 0*I ;

    % column 4 dF/dlDOP
    Jp{1,4} = -kappa_l*I ;
    Jp{2,4} = 0*I ;
    Jp{3,4} = 0*I;
    Jp{4,4} = TRdiv + kappa_l*I ;

    % right hand side of phosphate equations
    RHS = [kappa_g*DIPbar  ; ... 
           sparse(nwet,1)  ; ...
           sparse(nwet,1)  ; ...
           sparse(nwet,1)] ;

    %FD = mfactor(cell2mat(Jp)) ; 
    % Jacobian matrix
    J = cell2mat(Jp);
%%
    % Ftest = J*X - RHS;

    % % check that Ftest matches F
    % err = norm(F - Ftest);
    % if err > 1e-10
    %     error('Peqn: Jacobian check failed');
    % end
%%
end