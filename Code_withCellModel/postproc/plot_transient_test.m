

load('../output/PNAS2025_transient/transient_from_reoptNature_with_dop_GM15_npp1_CTL_He_PC.mat')

%% 
addpath('../../DATA/BGC_24layer/')
    addpath('../../DATA/BGC_2023Nature/')
    addpath('../utils/')

OperName = sprintf('OCIM2_CTL_He');
    load(OperName,'output') ;
    M3d = output.M3d;
    grd = output.grid;

    ilnd = find(M3d(:) == 0)    ;
iwet = find(M3d(:))         ;
nwet = length(iwet)         ;
dAt  = grd.DXT3d.*grd.DYT3d;
dVt  = dAt.*grd.DZT3d;

M3dsurf = M3d ;             % make surface mask (ocn grid cells in contact with atm) 
	M3dsurf(:,:,2:end) = 0 ;
	Msurf = M3dsurf(iwet);
	isrf = find(M3dsurf(iwet)) ;

spa = 60*60*24*365;
spd = 60*60*24;
%% check output

tindx = 2;

par.DIP = OUT.P(1:nwet,tindx);
par.POP = OUT.P(1*nwet+1:2*nwet,tindx);
par.DOP = OUT.P(2*nwet+1:3*nwet,tindx);
par.DOPl= OUT.P(3*nwet+1:4*nwet,tindx);


par.DIC    = OUT.C(1:nwet,tindx);
        par.POC    = OUT.C(1*nwet+1:2*nwet,tindx);
        par.DOC    = OUT.C(2*nwet+1:3*nwet, tindx);
        par.PIC    = OUT.C(3*nwet+1:4*nwet, tindx);
        par.ALK    = OUT.C(4*nwet+1:5*nwet, tindx);
        par.DOCl   = OUT.C(5*nwet+1:6*nwet, tindx);
        par.DOCr   = OUT.C(6*nwet+1:7*nwet, tindx);
        par.pco2atm= OUT.C(7*nwet+1, tindx);  

    DIC  = M3d+nan ;  DIC(iwet)  = par.DIC ;
    POC  = M3d+nan ;  POC(iwet)  = par.POC ;
    DOC  = M3d+nan ;  DOC(iwet)  = par.DOC ;
    PIC  = M3d+nan ;  PIC(iwet)  = par.PIC ;
    ALK  = M3d+nan ;  ALK(iwet)  = par.ALK ;
    DOCl = M3d+nan ;  DOCl(iwet) = par.DOCl ;
    DOCr = M3d+nan ;  DOCr(iwet) = par.DOCr ;
        % 

        sum(par.POP(isrf).*dVt(iwet(isrf)))

        %% 
 XX = load('../output/PNAS2025_transient/transient_test_6h_2d_5d_from_reoptNature_with_dop_GM15_npp1_CTL_He_PC_diagnostics.mat')
