% SetUp_v2.m
% uses input fields from DATA/BGC_2023Nature 
% but does not smooth po4 and tempobs data
% - load 2023Nature WOA files for po4 no3 temp sal
% - using old co2syspar for Fsea2air
global GC GO iter
iter = 0       ;
spd  = 24*60^2 ;
spa  = 365*spd ;
on   = true    ;
off  = false   ;
% addpath according to opterating system
% if GridVer_z ==48
%   addpath('../../DATA/BGC_48layer/')
% else
    addpath('../../DATA/BGC_24layer/')
    addpath('../../DATA/BGC_2023Nature/')
    addpath('../utils/')
% end

if GridVer == 91 
    switch(operator)
      case 'A'
        TRdivVer = 'CTL_He'             ;
      case 'B'
        TRdivVer = 'CTL_noHe'           ;
      case 'C'
        TRdivVer = 'KiHIGH_He'          ;
      case 'D'
        TRdivVer = 'KiHIGH_noHe'        ;
      case 'E'
        TRdivVer = 'KvHIGH_KiLOW_He'    ;
      case 'F'
        TRdivVer = 'KvHIGH_KiLOW_noHe'  ;
      case 'G'
        TRdivVer = 'KiLOW_He'           ;
      case 'H'
        TRdivVer = 'KiLOW_noHe'         ;
      case 'I'
        TRdivVer = 'KvHIGH_He'          ;
      case 'J'
        TRdivVer = 'KvHIGH_noHe'        ;
      case 'K'
        TRdivVer = 'KvHIGH_KiHIGH_noHe' ;
      case 'L'
        TRdivVer = 'CTL_He_48layer'     ;
    end 
end 

fprintf('Transport version: % s \n', TRdivVer)
if par.Cmodel == on
    fprintf('-------- C model is on -------- \n')
    fprintf('DOP scaling factor is %2.2e \n', par.dopscale)
    fprintf('DIP scaling factor is %2.2e \n', par.dipscale)
    fprintf('DIC scaling factor is %2.2e \n', par.dicscale)
    fprintf('DOC scaling factor is %2.2e \n', par.docscale)
    fprintf('ALK scaling factor is %2.2e \n', par.alkscale)
    fprintf('O2 scaling factor is %2.2e \n', par.o2scale)

  % set C2P function Type
  fprintf('------- Stoichiometry function ------------- \n')
  switch(par.C2Pfunctiontype)
    case 'P'
      par.C2P_PO4model = on;
      par.C2P_Tzmodel = off;
      par.Cellmodel = off;
      par.C2P_constant = off;
      par.C2P_loadprescribed = off;
    case 'C'
      par.Cellmodel = on;
      par.C2P_PO4model = off;
      par.C2P_Tzmodel = off;
      par.C2P_constant = off;
      par.C2P_loadprescribed = off;
    case 'T'
      par.C2P_Tzmodel = on;
      par.C2P_PO4model = off;
      par.Cellmodel = off;
      par.C2P_constant = off;
      par.C2P_loadprescribed = off;
    case 'R'
      par.C2P_constant = on;
      par.C2P_Tzmodel = off;
      par.C2P_PO4model = off;
      par.Cellmodel = off;
      par.C2P_loadprescribed = off;
    case 'L'
      par.C2P_constant = off;
      par.C2P_Tzmodel = off;
      par.C2P_PO4model = off;
      par.Cellmodel = off;
      par.C2P_loadprescribed = on;
  end

  if par.Cellmodel == on
    fprintf('-- Cell stoichiometry model is on  \n')
    if ~isfield(par,'dynamicP')
      fprintf('   -- default: Cell model depends on observed nutrient fields \n')
      par.dynamicP = off;
    elseif par.dynamicP == on
      fprintf('   -- Cell model depends on modelled DIP \n')
    else
      fprintf('   -- Cell model depends on observed nutrient fields \n')
    end
    % check C:P parameter flags
    if any([par.opt_cc, par.opt_dd, par.opt_ccT, par.opt_ddT])
      fprintf('Resetting opt_ccT, opt_ddT, opt_cc, and opt_dd to off ; cannot optimize linear C2P function parameters when cell model is on \n')
      par.opt_ccT = off;
      par.opt_ddT = off;
      par.opt_cc = off;
      par.opt_dd = off;
    end

  elseif par.C2P_Tzmodel == on
    fprintf('-- P:C is a linear function of WOA observed Temperature (normalized) \n')
    % check C:P parameter flags
    if any([par.opt_cc, par.opt_dd])
      fprintf('Resetting opt_cc and opt_dd to off ; cannot optimize cc and dd when C2P is a function of temperature \n')
      par.opt_cc = off;
      par.opt_dd = off;
    end

  elseif par.C2P_constant == on 
    fprintf('--- P:C is a constant value globally  \n')
    % check C:P parameter flags
    if any([par.opt_ccT, par.opt_ddT, par.opt_cc])
      fprintf('Resetting opt_cc, opt_ccT and opt_ddT to off ; only use opt_dd optimize a constant C2P value \n')
      par.opt_ccT = off;
      par.opt_ddT = off;
      par.opt_cc  = off;
    end

  elseif isfield(par,'C2P_loadprescribed') & par.C2P_loadprescribed == on
    fprintf('P:C pattern is prescribed. loaded from file: %s \n',par.fc2pload)
    % check C:P parameter flags
    if any([par.opt_cc, par.opt_dd, par.opt_ccT, par.opt_ddT])
      fprintf('Resetting opt_ccT, opt_ddT, opt_cc, and opt_dd to off ; cannot optimize linear C2P function parameters when using prescribed c2p \n')
      par.opt_ccT = off;
      par.opt_ddT = off;
      par.opt_cc = off;
      par.opt_dd = off;
    end
    % load C2P map
    % fprintf('C:P is optimal C:P from Sullivan et al 2024 (GBC) \n')
    % par.fc2pload = '../../DATA/BGC_24layer/C2Puptake_CellModel_opt_GBC2024.mat';
    load(par.fc2pload,'C2Puptake')
    par.C2Pprescribed = C2Puptake;
    par.p2c = 1./C2Puptake;

  else
    if ~isfield(par,'dynamicP')
      par.dynamicP = off;
      fprintf('-- P:C is a linear function of WOA observed phosphate  \n')
    elseif par.dynamicP == on
      fprintf('-- P:C is a linear function of modelled DIP  \n')
    else
      fprintf('-- P:C is a linear function of WOA observed phosphate \n')
    end
    % check C:P parameter flags
    if any([par.opt_ccT, par.opt_ddT])
      fprintf('Resetting opt_ccT and opt_ddT to off ; cannot optimize ccT and ddT when C2P is a function of phosphate \n')
      par.opt_ccT = off;
      par.opt_ddT = off;
    end
  end

else
  fprintf('--- Carbon cycle model is OFF ------ \n')
  par.C2P_constant = off;
  par.C2P_Tzmodel = off;
  par.C2P_PO4model = off;
  par.Cellmodel = off;
end  % end Cmodel
fprintf('\n')


if par.Omodel == on
    fprintf('-------- O model is on -------- \n')
    if isfield(par,'fo2cload')
      fprintf('Load O2C pattern from file: %s \n',par.fo2cload)
      load(par.fo2cload,'O2C')
      par.O2Cprescribed = O2C;
      fprintf('Resetting opt_O2C_T opt_rO2C to off ; cannot optimize O2C function parameters when using prescribed O2C \n')
      % tunable parameters;
      par.opt_O2C_T = off ;
      par.opt_rO2C  = off ;
    end
end 
if par.Simodel == on
    fprintf('-------- Si model is on -------- \n')
end 
fprintf('\n')

% ----- Load 24 layer Data -----
    OperName = sprintf('OCIM2_%s',TRdivVer);
    load(OperName,'output') ;
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
    %
    fname_npp = 'biopump_model_output_Nowicki.nc';
    fprintf('load NPP from %s ... \nNPP version = %d     (CbPM (1) or CAFE (2)) \n',fname_npp,par.nppVer)
    NPP = ncread(fname_npp,'NPP'); % 1 = CbPM; 2 = CAFE
    npp = NPP(:,:,par.nppVer); % mmol/m2/yr 
    %
    load M3d91x180x24.mat MSKS
    % WOA13 data
      load Sobs_91x180x24.mat Sobs    % Salinity obs for PME and Fsea2air
      salobs = Sobs ; clear Sobs;
      load po4obs_91x180x24.mat       % WOA PO4 climatological obs [units: umol/kg]
      DIP_obs = po4obs;  clear po4obs;
      %load no3obs_91x180x24.mat       % WOA NO3 clim obs [units: umol/kg]
      %DIN_obs = no3obs; clear no3obs;
      load tempobs_91x180x24.mat
      load Siobs_91x180x24.mat Siobs  % not needed for cell model
        Si_obs = Siobs; clear Siobs; 
    % WOA18 data
    % load TS_WOA_91x180x24.mat tempobs salobs % WOA temperature & salinity
    load O2_Nut_WOA_91x180x24.mat DIN_obs ; %O2_obs Si_obs DIN_obs DIP_obs % WOA O2 Si DIN DIP observations
        DIN_obs = DIN_obs*1.025; %[mmol/m^3]

    load PME_TS_91x180x24.mat pme % calculated from transport operator and salt fields (alternate to running pme.m) 
    load DICant_91x180x24.mat
    load GLODAPv2_91x180x24raw.mat alkraw po4raw o2raw sio4raw % GLODAP Nutrient units = [umol/kg]
    % load co2syspar_91x180x24.mat co2syspar % created in make_datafile_24layer/make_co2syspar.m
    load co2syspar91.mat co2syspar % file from 2023 Nature

    %load DOMobs_91x180x24.mat DOPobs   % old data file
    %load DOCobs_clean_91x180x24.mat    % old data file with refractory component and outlier removed
    load DOPobs_91x180x24.mat DOPobs
    load DOCobs_Feb2022_91x180x24.mat
    dopraw = DOPobs; 
    docraw = DOCobs;
    
    load initCO_91x180x24.mat data %initial guess for C&O
    load splco2_mod_monthly.mat             % monthly CO2 data fir transient run
    load GLODAPv2_DIC_remove_cant.mat gc12new
    dic_initial = nanmean(gc12new,4);

    load PARobs_processed_91x180x24.mat PARobs
    % PARobs = load('annual_PAR_91x180.mat'); %PAR data in units of [Einstein m-2 d-1] (units converted for cell model later in SetUp)
	  % PARobs.par = PARobs.PAR;
	  % load Kd490_MODIS_91x180.mat		% diffuse attenuation [m^-1]

   
%---------------------- constants -------------------
par.spd = spd ;
par.spa = spa ;
[par.kw, par.P] = kw(M3d, grd) ;
par.rho   = 1025         ; % seawater density; pme파일과 동일하게
permil    = par.rho*1e-3 ; % from umol/kg to mmol/m3; 
par.permil = permil      ;
% transiant CO2 concentraion ;
par.year      = splco2_mod(:,1) ;     
par.pco2_air  = splco2_mod(:,2) ;     
par.pco2atm   = splco2_mod(1,2) ;     % Year:1765 pCO2:277.9541 --> needed transient run
par.co2syspar = co2syspar       ;

%------------------- model grid info ----------------
ilnd = find(M3d(:) == 0)    ;
iwet = find(M3d(:))         ;
nwet = length(iwet)         ;
dAt  = grd.DXT3d.*grd.DYT3d;
dVt  = dAt.*grd.DZT3d;
ARC  = MSKS.ARC     ;
MED  = MSKS.MED     ;
PAC  = MSKS.PAC     ;
ATL  = MSKS.ATL     ;
IND  = MSKS.IND     ;
iarc = find(ARC(:))         ;
imed = find(MED(:))         ;
ipac = find(PAC(:))         ;
iatl = find(ATL(:))         ;
iind = find(IND(:))         ;
%
par.dAt   = dAt     ;
par.dVt   = dVt     ;
par.M3d   = M3d     ;
par.iwet  = iwet    ;
par.nwet  = nwet    ;
par.TRdiv = -TR     ;             % Convergence ----> Divergence
par.grd   = grd     ;
par.MSKS  = MSKS    ;
par.I = speye(nwet) ;

%------------------------  data  -----------------------
dicraw = dic_initial;           % ---> Transient run 통해 다시 re-optimized 필요
% get rid of mediterranean observations 
%     Med is disconnected from the ocean in the model because of course model resolution 
        % 기존 weilei 모델에서는 arc med자료 뺐음. 이번에는 다 포함시켜 보자.
fprintf('Removing mediterranean and arctic data from docraw, dicraw, alkraw ...\n')
docraw(imed)  = nan ;
alkraw(imed)  = nan ; 
dicraw(imed)  = nan ;
%po4raw(imed)  = nan ; 
%o2raw(imed)   = nan ;
% keep arctic observations
docraw(iarc)  = nan ;    
alkraw(iarc)  = nan ;
dicraw(iarc)  = nan ;
% po4raw(iarc)  = nan ;
% o2raw(iarc)   = nan ;

% Remove outliers from DOC observations
%   GLODAP climatologies tend to be less noisy, so outlier removal might be unnecessary
%   However, since parameter optimization is highly sensitive to outliers,
%   excluding extremes helps achieve the best fit to the majority of the data,
%   without undue influence from outliers
%   Note: Since rmOutliers removes a percentage of data, it removes a lot more 
%   data points from the GLODAP than DOC datasets
%
%   idea: we may want to remove outliers from minimization, but then keep all data 
%   when assessing how well the model fits the obs

fprintf('Removing outliers from docraw, po4raw, dicraw, alkraw, o2raw ...\n')
docraw = rmOutliers(docraw, M3d) ;
po4raw = rmOutliers(po4raw, M3d) ;
dicraw = rmOutliers(dicraw, M3d) ;
alkraw = rmOutliers(alkraw, M3d) ;
o2raw  = rmOutliers(o2raw, M3d)  ;

tempobs(tempobs(:)<-2.0) = -2.0   ; % Reset extreme cold temps to a minimum temperature threshold of -2C (seawater freezing point) 
DIP_obs(DIP_obs(:)<0.05) = 0.05   ;
po4raw(po4raw(:)<0.05)   = nan    ; % Remove DIP data below detection limit; 

if isfield(par,'SetUp_options') & par.SetUp_options.smoothP == 1
  fprintf('Smooth WOA DIP... \n')
  for ji = 1:24
     p2d = DIP_obs(:,:,ji);
     DIP_obs(:,:,ji) = smoothit(grd,M3d,p2d,3,1e5);   % 이게 필요한 이유...?
  end                                      %  ----> NPP field에 영향을 주려나?
else
  fprintf('Do not smooth WOA DIP... \n')
end                                   %  ----> NPP field에 영향을 주려나?
 
par.Temp     = tempobs       ;
par.Salt     = salobs        ; % Sobs     ;
par.DSi      = Si_obs        ;
par.po4obs   = DIP_obs       ;
par.no3obs   = DIN_obs       ;
par.o2raw    = o2raw         ;
%par.o2obs    = O2_obs        ;     %WOA O2 where is this used?
par.po4raw   = po4raw        ;
par.DOCobs   = docraw        ;
par.dicraw   = dicraw        ; 
par.alkraw   = alkraw*permil ;
par.dopraw   = dopraw - 0.05 ;    % To remove refractory component of DOP, subtract the mean deep ocean (>1000m) [DOP] (=0.05 μM, in DOPv2021 database (Liang et al., 2022) and at ALOHA (Foreman et al., 2019); following method of Liang et al., 2023)) 
par.docraw   = docraw ;
par.PARobs = PARobs;
%-------------------- prepare virtual flux -------------------
% PME part;
% [modT,modS,pme] = pme(par) ; 
par.pme = pme; %pme_new(iwet) ;   
junk = M3d ;
junk(:,:,2:end) = 0 ;
isrf = find(junk(iwet)) ;
sdic = find(par.dicraw(iwet(isrf)) > 0);
salk = find(par.alkraw(iwet(isrf)) > 0);
% surface mean concentraions
par.sDICbar = sum(par.dicraw((iwet(isrf(sdic)))).* ...
                  dVt(iwet(isrf(sdic))))./sum(dVt(iwet(isrf(sdic))));
par.sALKbar = sum(par.alkraw((iwet(isrf(salk)))).* ...
                  dVt(iwet(isrf(salk))))./sum(dVt(iwet(isrf(salk))));

%-------------------- normalize temperature --------------------
if isfield(par,'SetUp_options') & par.SetUp_options.smoothT == 1
  fprintf('Smooth WOA Temp... \n')
  for ji = 1:24
     t2d = par.Temp(:,:,ji); 
     par.Temp(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);        % smoothit ---> inpaint_nan하면 안돼?
  end                                                        % 일단 빼고 해보기.
else
  fprintf('Do not smooth WOA Temp... \n')
end

vT = par.Temp(iwet) ;                                     
% add + 1.0 to prevent from getting infinit kP or kC 
Tz = (vT - min(vT) + 1.0)./(max(vT) - min(vT)) ;
par.Tz = Tz      ;                                         % eqOcycle에서 수정해야함.
par.vT = vT      ;
Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
par.aveT   = nanmean(Tz3d(:,:,1:3),3) ;                    % tsnanmean하고 큰 차이?
                                                             % gp에서 run할때는
                                                             % nanmean으로 수정
                                                             % 1-3번째 layer의 T평균? 

% aveT: average of normalized T in top 3 layers. used in buildPFD to define a temperature dependent powerlaw exponent. power law assumes Temp is vertically uniform. 
%       % why normalized temp?
% vT: vectorized temperature. used in eqPcycle, eqCcycle, eqOcycle: for Q10 exponent for kC and kP
% Tz: normalized temperature. used in in eqCcycle or C2P_Tzmodel & in eqOcycle for O2C
                                                         

%-------------------- correct WOA o2 concentration --------------------
% o2obs_c = M3d*0;
% O2_obs(iwet) = O2_obs(iwet).*44.661;        % convert unit form [ml/l] to [umol/l].
% o2obs_c(iwet) = O2_obs(iwet).*1.009-2.523; % o2 correction based on Bianchi et al.(2012) [umol/l] .
% ineg = find(o2obs_c(:)<0);                % find out negative values and set them to zero.
% o2obs_c(ineg) = 0.1 ;
% par.o2obs = o2obs_c ;

%---------------------- prepare for restoring -----------------------
% calculating global mean DIP, ALK, and DSi concentraions for
% restoring 
idip = find(par.po4raw(iwet) > 0) ;
ialk = find(par.alkraw(iwet) > 0) ;
%isil = find(par.sio4raw(iwet)>0) ;

par.DIPbar = sum(par.po4raw(iwet(idip)).*dVt(iwet(idip)))/sum(dVt(iwet(idip))) ;
par.ALKbar = sum(par.alkraw(iwet(ialk)).*dVt(iwet(ialk)))/sum(dVt(iwet(ialk))) ;
%par.DSibar = sum(par.sio4raw(iwet(isil)).*dVt(iwet(isil)))/sum(dVt(iwet(isil)));

%for DIC (only needed for runs with Atmosphere box)
idic = find(par.dicraw(iwet)>0) ;
par.DICbar = sum(par.dicraw(iwet(idic)).*dVt(iwet(idic)))/sum(dVt(iwet(idic))) ;

% %------------------ Prepare Light field the model --------------------
% MOVED to cleanPARobs.m
% % par.nl = 2;
% PARobs = PARobs.par;

% % fill in missing values along coastlines
% PARsurf = inpaint_nans(PARobs);
% SURF = M3d(:,:,1);
% ilnd = find(SURF(:) == 0);
% PARsurf(ilnd) = NaN;
% PARsurf(PARsurf<=0) = min(PARobs(PARobs>0));

% % convert PAR [Einstein m^-2 d^-1] into units of [umol photon m^-2 s^-1] for cell model
% PARsurf = PARsurf*10^6/spd; % PAR at surface
% clear PARobs

% % extrapolate light to bottom of euphotic zone
% par.kI = 0.04;   % Light attenuation coefficient in seawater [m^-1]
% % Alternate:
%   % kI = 0.09*CHL.^0.4 ; % [m^-1]
%       % CESM light attenuation: average K for par wavelengths, plus absorbtion due to water; as a function of chlorophyll
%       % CHL = [mg/m^3]

%   %	kI = Kd490 * 73/2; % sum(grd.dzt(1:nl)) 
%       % CbPM model uses Satellite diffuse attenuation coefficient
%       %	 median mixed layer light level = surface irradiance * exp (-k490 * MLD/2)
% 	PAR        = 0*M3d;
% 	for ii=1:par.nl % only needed in euphotic zone for cell growth
% 		PAR(:,:,ii) = PARsurf.*exp(-par.kI*grd.zt(ii)); %PAR at mid depth of grid box [ii]
% 	end
  
%   par.PARobs = PAR;
% 	clear PAR

%-------------------- prepare NPP for the model ----------------------
% NPP unit (Nowicki) = (mmolC/m^2/yr)
% remove this P:C unit conversion. a constant stoichiometric scaling is implicit in alpha
if par.C2P_loadprescribed 
  fprintf('Sat NPP p2c conversion: 1/(prescribed C2Puptake field loaded from file) \n')
  %par.p2c = 1./C2Puptake; %set previously
elseif ~isfield(par,'SetUp_options') | par.SetUp_options.NPPp2c_type == 0
  fprintf('Sat NPP p2c conversion: GM15 \n')
  par.p2c = 0.006 + 0.0069*DIP_obs ;         
elseif par.SetUp_options.NPPp2c_type == 1
  fprintf('Sat NPP p2c conversion: 1/117 \n')
  par.p2c = (1/117) * M3d ;                
elseif par.SetUp_options.NPPp2c_type == 2
  fprintf('Sat NPP p2c conversion: optimal Cell Model C:P from Sullivan et al 2024 (GBC) \n')
  load('../../DATA/BGC_24layer/C2Puptake_CellModel_opt_GBC2024.mat','C2Puptake')
  par.p2c = 1./C2Puptake;
elseif par.SetUp_options.NPPp2c_type == 3
  fprintf('Sat NPP p2c conversion: Nature 2023 optimal C:P parameters : ')
  %fxhatload = '/DFS-L/DATA/primeau/hojons1/Nature2023_BGC_reoptimized/src_Nature_parameter_Megan/MSK91/CTL_He_PCO_Gamma0_kl12h_O5_POC2DIC_GM15_Nowicki_npp1_aveTeu_diffSig_O2C_uniEta_DICrmAnthro_2L_Pnormal_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00_xhat.mat' 
  tmp = load('/DFS-L/DATA/primeau/hojons1/Nature2023_BGC_reoptimized/src_Nature_parameter_Megan/MSK91/CTL_He_PCO_Gamma0_kl12h_O5_POC2DIC_GM15_Nowicki_npp1_aveTeu_diffSig_O2C_uniEta_DICrmAnthro_2L_Pnormal_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00_xhat.mat','xhat');
  fprintf('cc = %.3e ; dd = %.3e \n', tmp.xhat.cc, tmp.xhat.dd)
  par.p2c = 1./(tmp.xhat.cc * DIP_obs + tmp.xhat.dd);

end

inan = find(isnan(npp(:)) | npp(:) < 0) ;
npp(inan)  = 0  ;

par.nl = 2 ;                              % nl이 뭐지? NPP 발생하는 layer를 의미하나?
par.Lambda = M3d*0 ;
par.nppMSK = M3d*0 ;

% create mask for upper and lower ocean.
UMSK = M3d * 0 ;
DMSK = M3d * 0 ;
par.npp = repmat(npp, [1,1,par.nl]) * 0;
[nx,ny,nz] = size(M3d) ;
for jj = 1 : nx
    for ii = 1 : ny
        tmp = squeeze(M3d(jj, ii, :)) ;
        idp = length(find(tmp) == 1) ; %is this right?

        % if ( idp <= 12 )
        if ( idp <= 25 )
            UMSK(jj,ii,1:2) = 1 ;
            DMSK(jj,ii,3:end) = 1;
            par.nppMSK(jj,ii,1) = 2 ;
            par.nppMSK(jj,ii,2) = 2 ;

            par.npp(jj,ii,1)  = (1/2) * npp(jj,ii) / grd.dzt(1) * par.p2c(jj,ii,1);   % unit: (mmolP/m^3/yr)
            par.npp(jj,ii,2)  = (1/2) * npp(jj,ii) / grd.dzt(2) * par.p2c(jj,ii,2);

            Pnpp(jj,ii,1)  = par.npp(jj,ii,1) / spa ;                 % unit: (mmolP/m^3/s)
            Pnpp(jj,ii,2)  = par.npp(jj,ii,2) / spa ;

            Cnpp(jj,ii,1)  = (1/2) * npp(jj,ii) / grd.dzt(1) / spa ;  % unit: (mmolC/m^3/s)
            Cnpp(jj,ii,2)  = (1/2) * npp(jj,ii) / grd.dzt(2) / spa ;
                        
            par.Lambda(jj,ii,1) = 1./(1e-6+DIP_obs(jj,ii,1)) ;      % unit: [1/(mmolP/m^3)]
            par.Lambda(jj,ii,2) = 1./(1e-6+DIP_obs(jj,ii,2)) ;
        else 
            UMSK(jj,ii,1:3) = 1 ;
            DMSK(jj,ii,4:end) = 1;
            par.nppMSK(jj,ii,1) = 3 ;
            par.nppMSK(jj,ii,2) = 3 ;
            par.nppMSK(jj,ii,3) = 3 ;

            par.npp(jj,ii,1)  = (1/3) .* npp(jj,ii) / grd.dzt(1) * par.p2c(jj,ii,1);
            par.npp(jj,ii,2)  = (1/3) .* npp(jj,ii) / grd.dzt(2) * par.p2c(jj,ii,2);
            par.npp(jj,ii,3)  = (1/3) .* npp(jj,ii) / grd.dzt(3) * par.p2c(jj,ii,3);

            Pnpp(jj,ii,1)  = par.npp(jj,ii,1) / spa ;
            Pnpp(jj,ii,2)  = par.npp(jj,ii,2) / spa ;
            Pnpp(jj,ii,3)  = par.npp(jj,ii,3) / spa ;

            Cnpp(jj,ii,1)  = (1/3) * npp(jj,ii) / grd.dzt(1) / spa ;
            Cnpp(jj,ii,2)  = (1/3) * npp(jj,ii) / grd.dzt(2) / spa ;
            Cnpp(jj,ii,3)  = (1/3) * npp(jj,ii) / grd.dzt(3) / spa ;
            
            par.Lambda(jj,ii,1) = 1./(1e-6+DIP_obs(jj,ii,1)) ; 
            par.Lambda(jj,ii,2) = 1./(1e-6+DIP_obs(jj,ii,2)) ;
            par.Lambda(jj,ii,3) = 1./(1e-6+DIP_obs(jj,ii,3)) ;
        end 
    end
end

[par.Pnpp, par.Cnpp] = deal(M3d*0) ;
par.Pnpp(:,:,1:par.nl) = Pnpp; % unused
par.Cnpp(:,:,1:par.nl) = Cnpp;

UMSK(ilnd) = nan ;
DMSK(ilnd) = nan ; 
par.UM = d0(UMSK(iwet)) ;
par.DM = d0(DMSK(iwet)) ;
par.WM = d0(M3d(iwet))  ;
%---------------------------- end ---------------------------------- 
