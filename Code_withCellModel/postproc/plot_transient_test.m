

%load('../output/PNAS2025_transient/transient_test_PC_Atm_from_reoptNature_with_dop_GM15_npp1_CTL_He_PC_diagnostics.mat')

Xsteady = load('../output/PNAS2025_transient/transient_test_steadystate_PC_noAtm_from_reoptNature_with_dop_GM15_npp1_CTL_He_PC_diagnostics.mat')
Xoif = load('../output/PNAS2025_transient/transient_test_OIF_PC_noAtm_from_reoptNature_with_dop_GM15_npp1_CTL_He_PC_diagnostics.mat')
XequalbkCP = load('../output/PNAS2025_transient/transient_test_OIF_equalbkCP_PC_noAtm_from_reoptNature_with_dop_GM15_npp1_CTL_He_PC_diagnostics.mat')


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

%% HNLC regions map (publication-quality)
% build masks (as before)
MSKS.SO = M3d; % initialize mask of all ocean grid points
MSKS.SO(find(grd.yt>-45),:,:) = 0;

MSKS.EqP = M3d;
MSKS.EqP(find(grd.yt>5),:,:) = 0;
MSKS.EqP(find(grd.yt<-5),:,:) = 0;
MSKS.EqP(:,find(grd.xt>270),:) = 0;
MSKS.EqP(:,find(grd.xt<140),:) = 0;

MSKS.NP = M3d;
MSKS.NP(find(grd.yt>65),:,:) = 0;
MSKS.NP(find(grd.yt<45),:,:) = 0;
MSKS.NP(:,find(grd.xt<140),:) = 0;
MSKS.NP(:,find(grd.xt>220),:) = 0;

MSKS.HNLC = (MSKS.SO | MSKS.EqP | MSKS.NP);

% prepare surface (lon/lat) and masks
%[Lon,Lat] = meshgrid(grd.xt, grd.yt);    % Lon,Lat are size (ny x nx)
landMask = (M3d(:,:,1) == 0);            % true = land
hnlcMask = (MSKS.HNLC(:,:,1) > 0);       % true = HNLC at surface

% create an RGB background image (ocean / land colors)
ny = numel(grd.yt); nx = numel(grd.xt);
oceanColor = [0.86 0.92 0.98];   % pale blue
landColor  = [0.94 0.89 0.82];   % pale beige
bg = repmat(reshape(oceanColor,1,1,3), ny, nx);
for k = 1:3
    ch = bg(:,:,k);
    ch(landMask) = landColor(k);
    bg(:,:,k) = ch;
end

% plot
fig = figure('Units','pixels','Position',[100 100 1200 600]);
ax = axes('Parent',fig);
% draw background as a texture-mapped surface (no mapping toolbox required)
hSurf = surface('XData',grd.XT,'YData',grd.YT,'ZData',zeros(size(grd.XT)),...
    'CData',bg,'FaceColor','texturemap','EdgeColor','none','Parent',ax);
hold(ax,'on');

% overlay HNLC regions as semi-transparent filled squares at cell centers
% use square markers to reflect grid cells
idx = find(hnlcMask);
if ~isempty(idx)
    sc = scatter(grd.XT(idx), grd.YT(idx), 80, [0.85 0.33 0.10], 's', 'filled', ...
        'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'none', 'Parent', ax);
else
    sc = plot(NaN,NaN); % placeholder for legend if no points
end

% draw coastline contour from landMask for crisp boundary
contour(grd.XT, grd.YT, double(landMask), [0.5 0.5], 'k', 'LineWidth', 0.8, 'Parent', ax);

% aesthetics
axis(ax,'equal');
xlim(ax,[min(grd.xt) max(grd.xt)]);
ylim(ax,[min(grd.yt) max(grd.yt)]);
set(ax,'YDir','normal');   % ensure lat increases upward
box(ax,'on');
xlabel(ax,'Longitude');
ylabel(ax,'Latitude');
set(ax,'FontSize',14, 'TickDir','out');

% title and legend
%regions (areas where production rate was increased)
title(ax,'Fertilization experiment: HNLC Map','FontSize',16,'FontWeight','normal');
if exist('sc','var') && ~isempty(sc) && isvalid(sc)
    legend(sc,'Fertilized regions','Location','southoutside','Orientation','horizontal','Box','off');
end

% tidy for publication: remove minor visual clutter
set(fig,'Color','w');
colormap(ax,[]);
hold(ax,'off');

% optionally save a high-resolution PNG (uncomment to use)
% print(fig,'HNLC_regions_map','-dpng','-r300');


%%
Times = Diags.Tout/(60*60*24*365);
indx = find(Times<=79);
figure; plot(Times(indx),Diags.pco2atm(indx)); ylabel('pCO2atm')

%% total DIP
tmaxyr = 200;
indx = find(Times<=tmaxyr);
figure; plot(Times(indx),Diags.totalDIP(indx),'LineWidth',2); ylabel('totalDIP'); xlabel('Time (years)')
xlim([0,tmaxyr])

%%
indx = 1:length(Times);
figure; plot(Times(indx),Diags.totalDIPsurf(indx),'Linewidth',2); ylabel('totalDIP euphotic'); xlabel('Time (years)')

%% totalDIC
tmaxyr = 190;
indx = find(Times<=tmaxyr);
figure; plot(Times(indx),Diags.totalDIC(indx),'LineWidth',2); ylabel('totalDIC'); xlabel('Time (years)')
xlim([0,tmaxyr])

%% PNPP
indx = 1:length(Times); %find(Times<=10);
figure; 
plt1 = plot(Times(indx),Diags.PNPP(indx),'linewidth',2); ylabel('PNPP (Pg P/yr)'); xlabel('Time (years)'); title('NPP (Pg P/yr)','FontSize',14)
%% CNPP
Diags = Xoif.Diags;
indx = find(Times<=10);
figure; plot(Times(indx),Diags.CNPP_nolabile(indx),'linewidth',2,'DisplayName','1y HNLC fertilized x1.5'); ylabel('non-labile NPP (Pg C/yr)'); xlabel('Time (years)'); 
title('NPP excluding labile DOC (Pg C/yr)','FontSize',14)
hold on;
plot(Times(indx),ones(size(indx))*28.836,'--k','linewidth',2,'DisplayName','Steady-State')
% steady state: CNPP_nolab=28.836 PgC/yr
lgd = legend;
lgd.FontSize = 12;


hold on
plot(Times(indx),XequalbkCP.Diags.CNPP_nolabile(indx),'linewidth',2,'DisplayName','equal C &P remin; 1y HNLC fertilized x1.5'); ylabel('non-labile NPP (Pg C/yr)'); xlabel('Time (years)'); 

plot(Times(indx),Xsteady.Diags.CNPP_nolabile(indx),'linewidth',2,'DisplayName','steady state'); ylabel('non-labile NPP (Pg C/yr)'); xlabel('Time (years)'); 

%% CNPP

Diags = XequalbkCP.Diags;

% [0.7 0.7 0.7]
colors.medblue = [0 0.4470 0.7410];
%indx = find(Times<=10);
indx = 1:length(Times);
figure; 
% steady state: CNPP_nolab=28.836 PgC/yr
plot(Times(indx),ones(size(indx))*28.836,'-','Color','k','linewidth',2,'DisplayName','Steady-State')
hold on;
plot(Times(indx),Diags.CNPP_nolabile(indx),'-','Color',colors.medblue,'linewidth',2,'DisplayName','1y HNLC fertilized x1.5'); ylabel('non-labile NPP (Pg C/yr)'); xlabel('Time (years)'); 
%plot(Times(indx),Diags.CNPP_nolabile(indx),'.b','MarkerFaceColor','b','linewidth',2,'DisplayName','1y HNLC fertilized x1.5','MarkerSize',10); ylabel('non-labile NPP (Pg C/yr)'); xlabel('Time (years)'); 
title('NPP excluding labile DOC (Pg C/yr)','FontSize',14)
axis tight;

lgd = legend;
lgd.FontSize = 12;

%% plot the difference between Diags.CNPP_nolabile and steady state CNPP
indx = find(Times>1);
figure;
pltdata = Diags.CNPP_nolabile(indx) - 28.836;
plot(Times(indx),pltdata,'-','Color',colors.medblue,'linewidth',2,'DisplayName','1y HNLC fertilized x1.5'); ylabel('non-labile NPP Difference from steady state (Pg C/yr)'); xlabel('Time (years)'); 
title('Transient change in NPP from steady-state excluding labile DOC (Pg C/yr)','FontSize',14)
%axis tight
ylim([-0.55, 0.05])
xlim([0,200])

%% plot total C export
Times = Xoif.Diags.Tout/spa;
indx = find(Times<=10);
figure;

pltdata = Xoif.Diags.totalCexport_darkremin(indx);
plot(Times(indx),pltdata,'-','Color',colors.medblue,'linewidth',2,'DisplayName','optimized model'); 
ylabel('C export (Pg C/yr)'); xlabel('Time (years)'); 
title('C export: 1y HNLC fertilized x1.5')
hold on;

pltdata = XequalbkCP.Diags.totalCexport_darkremin(indx);
plot(Times(indx),pltdata,'-','Color','m','linewidth',2,'DisplayName','equal P & C remin'); 

pltdata = Xsteady.Diags.totalCexport_darkremin(indx);
plot(Times(indx),pltdata,'-','Color','k','linewidth',2,'DisplayName','steady state optimized model'); 

legend

%% plot total C export 200 years
Times = Xoif.Diags.Tout/spa;
indx = find(Times<=202);
figure; hold on;

pltdata = Xsteady.Diags.totalCexport_darkremin(indx);
plot(Times(indx),pltdata,'-','Color','k','linewidth',2,'DisplayName','steady state optimized model'); 

pltdata = Xoif.Diags.totalCexport_darkremin(indx);
plot(Times(indx),pltdata,'-','Color',colors.medblue,'linewidth',2,'DisplayName','1y HNLC fertilized x1.5'); 
ylabel('C export (Pg C/yr)'); xlabel('Time (years)'); 
title('C export: 1y HNLC fertilized x1.5')


%pltdata = XequalbkCP.Diags.totalCexport_darkremin(indx);
%plot(Times(indx),pltdata,'-','Color','m','linewidth',2,'DisplayName','equal P & C remin'); 


legend
xlim([0, 200])

%% plot delta C export from steady state

Times = Xoif.Diags.Tout/spa;
tmaxyr = 10;
indx = find(Times<=tmaxyr);
figure; hold on;

pltdata = Xoif.Diags.totalCexport_darkremin(indx) - Xsteady.Diags.totalCexport_darkremin(indx);
plot(Times(indx),zeros(size(pltdata)),'k-','LineWidth',2,'DisplayName','Zero')
plot(Times(indx),pltdata,'-','Color',colors.medblue,'linewidth',2,'DisplayName','fertilized'); 
ylabel('Change in C export (Pg C/yr)'); xlabel('Time (years)'); 
title('Change in organic C export: 1y HNLC fertilized x1.5')
%legend
xlim([0, tmaxyr])
grid on
%%
% at what time, does CNPP_nolabile reach CNPP_steadystate?

%% calculate difference between fertilized and no fertilized run
% net integrated change in export production

%% from Diags, check that total P content of ocean stays the same
Natm    = 1.773e20  ;
totalPocn = Diags.totalDIP + Diags.totalPOP + Diags.totalDOP + Diags.totalDOPl; 
totalCocn = Diags.totalDIC + Diags.totalPOC + Diags.totalDOC + Diags.totalDOCl + Diags.totalDOCr + Diags.totalPIC;
totalCO2atm = Diags.pco2atm * Natm * 1e-6 * 12 * 1e-15; % Pg C
totalCsystem = totalCocn + totalCO2atm; 

fprintf('Integrated total C range: %10.3e Pg C  -  %10.3e Pg C \n',min(totalCsystem),max(totalCsystem));
fprintf('Integrated total P range: %10.3e Pg P  -  %10.3e Pg P \n',min(totalPocn),max(totalPocn));

%% make a plot of totalCsystem
indx = find(Times<=90);
figure; plot(Times(indx),totalCsystem(indx)); ylabel('total C system (Pg C)'); xlabel('Time (years)')  
%%
indx = 1:length(Tout);
figure; plot(Times(indx),totalPocn(indx),'b-','linewidth',2); ylabel('total P ocean (Pg P)'); xlabel('Time (years)')
title('Total P inventory in ocean', 'FontSize',14)

%% 
% add a diagnostic for C export? or POC export specifically?

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
