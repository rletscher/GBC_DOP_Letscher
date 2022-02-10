
% 30 free parameters

% prescribe some constants.
dpa = 365;      % days per year;
spd = 24*60^2;  % seconds per day;
spa = dpa*spd;  % seconds per year;

% OCIM transport operator including grid (grid) and mask (M3d) and
% advection and diffusion transport operator (TR);   
load OCIM2/OCIM2_CTL_He.mat
TR = output.TR; % TR is the flux ***convergence***
TR = TR./spa;   % convert from yr-1 to s-1 units
load data/ao.mat % load grid-describing strcture ao
% 'M3d' is land-sea mask
M3d = output.M3d;
% 'grid' contains grid variables and metrics
grid = output.grid;
%'msk' is a stucture with fields representing the indices of points
msk =output.msk;
% 'icon' are tracer grid points 
iocn = msk.pkeep;
% calculate volume in each grid 
VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d; % m^3


%%%%%% create omega - 9 region mask %%%%%%
i=length(iocn);
load masks_DOP.mat
omega(:,1)=region1(iocn); omega(:,2)=region2(iocn); omega(:,3)=region3(iocn);
omega(:,4)=region4(iocn); omega(:,5)=region5(iocn); omega(:,6)=region6(iocn);
omega(:,7)=region7(iocn); omega(:,8)=region8(iocn); omega(:,9)=region9(iocn);

%%%%%% load global PO4 concentration field %%%%%%
load data/WOA13/WOAPO4.mat % load WOA13 PO4 observations
po4obs = WOAPO4;

%%%%%% load DOM observations %%%%%%
dataset = 'DOPv2021v4.csv';
d=4; % set for case: 4=DOP 
DOM_obs = get_bottle_dom(grid,M3d,d,dataset); %mmol m^-3
DOM_obs = DOM_obs(:,:,:)-0.05; % subtract out DOPr 0.05 µM
% eliminate DOM observtion < 0
DOM_obs(DOM_obs<0)=nan;
% refill land with nan
DOM_obs(M3d(:)==0)=nan;



grd   = grid;         % grid info. for OCIM model;
iwet  = find(M3d(:)); % find index of wet grid box;
nwet  = length(iwet); % number of wet grid boxes;
I     = speye(nwet);  % build an identity matrix;

% OCIM model info. loaded from ao;
parm.M3d = M3d;   % land ocean mask;
parm.TRdiv = -TR; % advection/diffusion transport operator [s^-1];
parm.grd = grd;   % OCIM grid info
dVt = grid.DZT3d; % OCIM grid box height
parm.dVt = dVt;   % depth level thickness [m]


% wrap data into structure for easy passing to other functions;
parm.p2c      = 0.006+0.0069*po4obs; % P:C ratio (GM15) used to scale NPP
parm.DIPbar   = sum(po4obs(iwet).*dVt(iwet))./sum(dVt(iwet));
parm.DIPstd   = std(po4obs(iwet));
parm.po4obs   = po4obs;
parm.DOPobs   = DOM_obs;
parm.omega    = omega;

ipo4 = find(~isnan(po4obs(iwet)));  % index for valid PO4 measurements
idop = find(~isnan(DOM_obs(iwet))); % index for valid DOP measurements

n_po4 = length(ipo4); % number of valid PO4 obs
n_dop = length(idop); % number of valid DOP obs


%%%%%%% prepare NPP for the model %%%%%%%%
%%%%%% load global NPP field %%%%%%%
load SeaWiFS_annual_CbPM_NPP_gC_m2_yr_onOCIM2.mat  % unit: g C/m2/yr
NPP = NPP(:,:).*(1000/12); % unit: mmol C/m2/yr
inan = find(isnan(NPP(:)));
NPP(inan) = 0;
parm.npp = NPP(:,:,:)./spa;
parm.Lambda = M3d*0;
% divide up NPP evenly in top 2 layers; convert C to P using GM15; convert
% m-2 to m-3; units = 1/meter*carbon
parm.Lambda(:,:,1) = 0.5*(1/grd.dzt(1))*parm.p2c(:,:,1)./(1e-9+po4obs(:,:,1));
parm.Lambda(:,:,2) = 0.5*(1/grd.dzt(2))*parm.p2c(:,:,2)./(1e-9+po4obs(:,:,2));
parm.Lambda(:,:,3:end) = 0;
%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% prescribed parameters %%%%%%%%%%%%%%
parm.kappa_g  = 1/(1e6*spa);  % geological restoring;
parm.kappa_p  = 1/(30*spd);   % POP remin to DIP with 30 day lifetime; Wang et al 2019
parm.v        = 1/(20*spa);   % min autotroph DOP uptake lifetime of 20 yr
%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% parameters need to be optimized %%%%%%%
sigma0    = 0.1;          % fraction of NPP directed to DOPsl; global model posterior
kappa0    = 1/(3.5*spa);  % remin rate of DOPsl converted to s-1; 2x faster than DOC Letscher et al 2015
km0      = 0.2;           % half sat constant for DOPsl autotrophic uptake µM; Mather et al 2008
alpha    = 2.0e-4;        % Pi uptake from NPP fitting parameter; constant s-1; global model posterior
beta     = 0.036;         % Pi uptake from NPP fitting parameter; exponent; global model posterior
bp       = 1.16;          % Martin's 'b' exponent for POP attenuation; global model posterior
%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% build regional parameter vector %%%%%%%%
sigma=zeros(9,1);
sigma(:)=sigma0;
kappa=zeros(9,1);
kappa(:)=kappa0;
km=zeros(9,1);
km(:)=km0;

% log transform parameters that have to be positive;
x0 = log([sigma; kappa; km; alpha; beta; bp;]);

nip = length(x0); % number of optimized parameters;

% order of second derivatives
in = 1:nip; % parameter indices
nin = length(in);
n = (nin+1)*nin/2;       % number of unique 2nd derivatives
                         
% the pos matrix is an nip x nip array with numbers from
% 1 to n in the upper triangular part. These numbers are used as
% the column index in a big array whos columns store the 2nd
% derivative of the model state with respect to the parameters
pos = tril(ones(nin),0); 
pos(pos==1) = 1:n;
pos = pos'; 
parm.in = in;
parm.nin = nin;
parm.pos = pos;


% pass the initial iterate and parameters to the optimizer
myfun = @(x) neglogpost_DOP(x,parm);
options = optimoptions(@fminunc,...
                       'Algorithm','trust-region',...
                       'GradObj','on',...
                       'Hessian','on',...
                       'Display','iter',...
                       'MaxFunEvals',6000, ...
                       'MaxIter',6000,...
                       'TolX',1e-9,...
                       'TolFun',1e-9,...
                       'DerivativeCheck','off',...
                       'FinDiffType','central',...
                       'PrecondBandWidth',Inf);

[xhat,fval,exitflag] = fminunc(myfun,x0,options);

[f,fx,fxx] = neglogpost_DOP(xhat,parm);

save('cost','f','fx','fxx','xhat')

    %%% cost.mat holds
    %%% f    = min value of cost function (scalar)
    %%% fx   = first deriv of cost func w.r.t. params (1x30)
    %%% fxx  = second deriv of cost func w.r.t. params (30x30)
    %%% xhat = log of optimized param values (30x1)

        



%redo calculation
p        = [];
p.sigma1 = exp(xhat(1));
p.sigma2 = exp(xhat(2));
p.sigma3 = exp(xhat(3));
p.sigma4 = exp(xhat(4));
p.sigma5 = exp(xhat(5));
p.sigma6 = exp(xhat(6));
p.sigma7 = exp(xhat(7));
p.sigma8 = exp(xhat(8));
p.sigma9 = exp(xhat(9));
p.kappa1 = exp(xhat(10));
p.kappa1 = 1/(p.kappa1*spa);
p.kappa2 = exp(xhat(11));
p.kappa2 = 1/(p.kappa2*spa);
p.kappa3 = exp(xhat(12));
p.kappa3 = 1/(p.kappa3*spa);
p.kappa4 = exp(xhat(13));
p.kappa4 = 1/(p.kappa4*spa);
p.kappa5 = exp(xhat(14));
p.kappa5 = 1/(p.kappa5*spa);
p.kappa6 = exp(xhat(15));
p.kappa6 = 1/(p.kappa6*spa);
p.kappa7 = exp(xhat(16));
p.kappa7 = 1/(p.kappa7*spa);
p.kappa8 = exp(xhat(17));
p.kappa8 = 1/(p.kappa8*spa);
p.kappa9 = exp(xhat(18));
p.kappa9 = 1/(p.kappa9*spa);
p.km1    = exp(xhat(19));
p.km2    = exp(xhat(20));
p.km3    = exp(xhat(21));
p.km4    = exp(xhat(22));
p.km5    = exp(xhat(23));
p.km6    = exp(xhat(24));
p.km7    = exp(xhat(25));
p.km8    = exp(xhat(26));
p.km9    = exp(xhat(27));
p.alpha  = exp(xhat(28));
p.beta   = exp(xhat(29));
p.bp     = exp(xhat(30));

fprintf('p.sigma1, p.sigma2, p.sigma3, p.sigma4, p.sigma5, p.sigma6, p.sigma7, p.sigma8, p.sigma9, p.kappa1, p.kappa2, p.kappa3, p.kappa4, p.kappa5, p.kappa6, p.kappa7, p.kappa8, p.kappa9, p.km1, p.km2, p.km3, p.km4, p.km5, p.km6, p.km7, p.km8, p.km9, p.alpha, p.beta, p.bp, %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f \n\n', p.sigma1, p.sigma2, p.sigma3, p.sigma4, p.sigma5, p.sigma6, p.sigma7, p.sigma8, p.sigma9, p.kappa1, p.kappa2, p.kappa3, p.kappa4, p.kappa5, p.kappa6, p.kappa7, p.kappa8, p.kappa9, p.km1, p.km2, p.km3, p.km4, p.km5, p.km6, p.km7, p.km8, p.km9, p.alpha, p.beta, p.bp);
  save p.mat p

    %%% p.mat holds
    %%% p    = optimized param values (30x1)



