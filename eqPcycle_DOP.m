function [P,Px,Pxx,parm] = eqPcycle_DOP(parm,ip,x)
% ip is the mapping from x to parameter names (see switch below)
% output: P is model prediction of DIP, POP, and DOP
% output: Px partial derivative of P model w.r.t. model parameters x
% output: Pxx hessian matrix of P model w.r.t. model parameters x
    
% unpack the parameters to be optimized
    if (nargin>1)
        for ik1 = 1:length(ip)
            switch ip(ik1)
              case 1
                parm.sigma1 = exp(x(ik1)); % fraction of NPP directed to DOPsl in region 1
              case 2
                parm.sigma2 = exp(x(ik1)); % fraction of NPP directed to DOPsl in region 2
              case 3
                parm.sigma3 = exp(x(ik1)); % fraction of NPP directed to DOPsl in region 3
              case 4
                parm.sigma4 = exp(x(ik1)); % fraction of NPP directed to DOPsl in region 4
              case 5
                parm.sigma5 = exp(x(ik1)); % fraction of NPP directed to DOPsl in region 5
              case 6
                parm.sigma6 = exp(x(ik1)); % fraction of NPP directed to DOPsl in region 6
              case 7
                parm.sigma7 = exp(x(ik1)); % fraction of NPP directed to DOPsl in region 7
              case 8
                parm.sigma8 = exp(x(ik1)); % fraction of NPP directed to DOPsl in region 8
              case 9
                parm.sigma9 = exp(x(ik1)); % fraction of NPP directed to DOPsl in region 9
              case 10
                parm.kappa1 = exp(x(ik1)); % DOPsl heterotrophic remin rate in region 1
              case 11
                parm.kappa2 = exp(x(ik1)); % DOPsl heterotrophic remin rate in region 2
              case 12
                parm.kappa3 = exp(x(ik1)); % DOPsl heterotrophic remin rate in region 3
              case 13
                parm.kappa4 = exp(x(ik1)); % DOPsl heterotrophic remin rate in region 4
              case 14
                parm.kappa5 = exp(x(ik1)); % DOPsl heterotrophic remin rate in region 5
              case 15
                parm.kappa6 = exp(x(ik1)); % DOPsl heterotrophic remin rate in region 6
              case 16
                parm.kappa7 = exp(x(ik1)); % DOPsl heterotrophic remin rate in region 7
              case 17
                parm.kappa8 = exp(x(ik1)); % DOPsl heterotrophic remin rate in region 8
              case 18
                parm.kappa9 = exp(x(ik1)); % DOPsl heterotrophic remin rate in region 9
              case 19
                parm.km1 = exp(x(ik1)); % DOPsl autrophic uptake half sat constant in region 1
              case 20
                parm.km2 = exp(x(ik1)); % DOPsl autrophic uptake half sat constant in region 2
              case 21
                parm.km3 = exp(x(ik1)); % DOPsl autrophic uptake half sat constant in region 3
              case 22
                parm.km4 = exp(x(ik1)); % DOPsl autrophic uptake half sat constant in region 4
              case 23
                parm.km5 = exp(x(ik1)); % DOPsl autrophic uptake half sat constant in region 5
              case 24
                parm.km6 = exp(x(ik1)); % DOPsl autrophic uptake half sat constant in region 6
              case 25
                parm.km7 = exp(x(ik1)); % DOPsl autrophic uptake half sat constant in region 7
              case 26
                parm.km8 = exp(x(ik1)); % DOPsl autrophic uptake half sat constant in region 8
              case 27
                parm.km9 = exp(x(ik1)); % DOPsl autrophic uptake half sat constant in region 9
              case 28
                parm.alpha = exp(x(ik1)); % npp scaling factor for DIP uptake rate
              case 29
                parm.beta = exp(x(ik1));  % npp scaling exponent for DIP uptake rate
              case 30
                parm.bp = exp(x(ik1));    % Martin exponent for POP remin
            end
        end
    end
   
    
  % regional parameters
    omega = parm.omega;              % 9 region masks in size iwet x 9
    
    % unpack some useful stuff
    M3d   = parm.M3d;
    TRdiv = parm.TRdiv;
    grd   = parm.grd;
    
    iwet  = find(M3d(:));         % wet point indices;
    nwet = length(iwet);        % number of wet points;
    
    I = speye(nwet);            % make an identity matrix;
    d0 = @(x) spdiags([x(:)],[0],length(x(:)),length(x(:)));

    
sigma1 = parm.sigma1;  
sigma2 = parm.sigma2; 
sigma3 = parm.sigma3; 
sigma4 = parm.sigma4; 
sigma5 = parm.sigma5; 
sigma6 = parm.sigma6; 
sigma7 = parm.sigma7; 
sigma8 = parm.sigma8; 
sigma9 = parm.sigma9;
% build 9 region sigvec (sigma vector) and sigmat (sigma matrix)
ss1=sigma1*omega(:,1); ss2=sigma2*omega(:,2); ss3=sigma3*omega(:,3);
ss4=sigma4*omega(:,4); ss5=sigma5*omega(:,5); ss6=sigma6*omega(:,6);
ss7=sigma7*omega(:,7); ss8=sigma8*omega(:,8); ss9=sigma9*omega(:,9);
sigvec = ss1+ss2+ss3+ss4+ss5+ss6+ss7+ss8+ss9;

sigmat = nan*M3d;
sigmat(iwet) = sigvec(:);
sigomega = d0(sigmat(iwet));
sigmat1=-1+sigmat;  % perform -1 addition before creating iwet x iwet sparse matrix 
sigomega1=d0(sigmat1(iwet));

kappa1 = parm.kappa1;
kappa2 = parm.kappa2;
kappa3 = parm.kappa3;
kappa4 = parm.kappa4;
kappa5 = parm.kappa5;
kappa6 = parm.kappa6;
kappa7 = parm.kappa7;
kappa8 = parm.kappa8;
kappa9 = parm.kappa9;
% build 9 region kapvec (kappa vector) and kapmat (kappa matrix)
kk1=kappa1*omega(:,1); kk2=kappa2*omega(:,2); kk3=kappa3*omega(:,3);
kk4=kappa4*omega(:,4); kk5=kappa5*omega(:,5); kk6=kappa6*omega(:,6);
kk7=kappa7*omega(:,7); kk8=kappa8*omega(:,8); kk9=kappa9*omega(:,9);
kapvec = kk1+kk2+kk3+kk4+kk5+kk6+kk7+kk8+kk9; 
kapmat = nan*M3d;
kapmat(iwet) = kapvec(:);
kapomega = d0(kapmat(iwet));

km1    = parm.km1;
km2    = parm.km2;
km3    = parm.km3;
km4    = parm.km4;
km5    = parm.km5;
km6    = parm.km6;
km7    = parm.km7;
km8    = parm.km8;
km9    = parm.km9;
% build 9 region kmvec (kmm)
kmm1=km1*omega(:,1); kmm2=km2*omega(:,2); kmm3=km3*omega(:,3);
kmm4=km4*omega(:,4); kmm5=km5*omega(:,5); kmm6=km6*omega(:,6);
kmm7=km7*omega(:,7); kmm8=km8*omega(:,8); kmm9=km9*omega(:,9);


alpha  = parm.alpha;       
beta   = parm.beta; 
bp     = parm.bp;



    % fixed parameters
    DIPbar  = M3d(iwet)*parm.DIPbar; % gobal arerage PO4 conc.[mmol m^-3]; 
    kappa_g = parm.kappa_g;          % PO4 geological restore const.[s^-1];
    kappa_p = parm.kappa_p;          % POP solubilization rate constant
    v       = parm.v;                % min autotroph DOP uptake lifetime of 20 yr;
    npp     = parm.npp;              % net primary production in mmol C m-2 s-1

    

    
    % build part of the biological DIP uptake operator
    Lambda = parm.Lambda;
    LAM = 0*M3d;
    LAM(:,:,1) = (npp.^beta).*Lambda(:,:,1);
    LAM(:,:,2) = (npp.^beta).*Lambda(:,:,2);
    L = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
    parm.L = L; 
    
    % build autotrophic uptake of DOP operator mu
    po4_tmp = parm.po4obs;        % WOA PO4 µM
    po4     = po4_tmp(iwet);            

    mu = d0( (v.*(kmm1+po4))./po4 + (v.*(kmm2+po4))./po4 + ...
             (v.*(kmm3+po4))./po4 + (v.*(kmm4+po4))./po4 + ...
             (v.*(kmm5+po4))./po4 + (v.*(kmm6+po4))./po4 + ...
             (v.*(kmm7+po4))./po4 + (v.*(kmm8+po4))./po4 + ...
             (v.*(kmm9+po4))./po4 ) ; 
    
    Pi = po4;               % need po4 vector for km gradient
    
    % build the sinking particulate flux divergence operator
    if (~isempty(find(ip==30)))
        % do this if b needs to be optimized
        [PFdiv,dPFDdb,d2PFDdb2] ...
            = buildPFD(M3d,bp,kappa_p,grd);
    else
        % do this if the  parameters don't need to be optimized
        PFdiv = buildPFD(M3d,bp,kappa_p,grd);  % particle flux divergence [s^-1];
    end
                   
    % build Jacobian matrix
    %    [                 dF/dDIP,          dF/dPOP,                     dF/dDOP]            
    Fp = [[  TRdiv+alpha*L+kappa_g*I,        -kappa_p*I,                  -kapomega.*I]; ...   % DIP eqn
          [        sigomega1*alpha*L,   PFdiv+kappa_p*I,                           -mu]; ...   % POP eqn            
          [        -sigomega*alpha*L,               0*I,        (TRdiv+kapomega.*I)+mu]];      % DOP eqn
    
    % right hand side of phosphate equations 
    Z=zeros(nwet,1);
    RHS = [ kappa_g*DIPbar;...
                         Z;...
                         Z];
    
    % dP/dt + Fp*P = RHS    
    % factoring Jacobian matrix
    FFp = mfactor(Fp); 
    % solve for P-cycle model state
    P = mfactor(FFp,RHS);    % the beauty of linear algebra!
    if (nargout>1)
        %
        %
        % Compute the gradient of the solution wrt the parameters
        %
        %
        DIP = P(1:nwet);
        POP = P(nwet+1:2*nwet);
        DOP = P(2*nwet+1:end);
        Fx = zeros(3*nwet,length(ip));
        for ik1 = 1:length(ip)
            switch ip(ik1)

              case 1 % sigma1
                Fx(:,ik1) =  exp(x(ik1))* ...
                    [Z;...
                     alpha*d0(omega(:,1))*(L*DIP);...
                     -alpha*d0(omega(:,1))*(L*DIP)];
              case 2 % sigma2
                Fx(:,ik1) =  exp(x(ik1))* ...
                    [Z;...
                     alpha*d0(omega(:,2))*(L*DIP);...
                     -alpha*d0(omega(:,2))*(L*DIP)];
              case 3 % sigma3
                Fx(:,ik1) =  exp(x(ik1))* ...
                    [Z;...
                     alpha*d0(omega(:,3))*(L*DIP);...
                     -alpha*d0(omega(:,3))*(L*DIP)];
              case 4 % sigma4
                Fx(:,ik1) =  exp(x(ik1))* ...
                    [Z;...
                     alpha*d0(omega(:,4))*(L*DIP);...
                     -alpha*d0(omega(:,4))*(L*DIP)];
              case 5 % sigma5
                Fx(:,ik1) =  exp(x(ik1))* ...
                    [Z;...
                     alpha*d0(omega(:,5))*(L*DIP);...
                     -alpha*d0(omega(:,5))*(L*DIP)];
              case 6 % sigma6
                Fx(:,ik1) =  exp(x(ik1))* ...
                    [Z;...
                     alpha*d0(omega(:,6))*(L*DIP);...
                     -alpha*d0(omega(:,6))*(L*DIP)];
              case 7 % sigma7
                Fx(:,ik1) =  exp(x(ik1))* ...
                    [Z;...
                     alpha*d0(omega(:,7))*(L*DIP);...
                     -alpha*d0(omega(:,7))*(L*DIP)];
              case 8 % sigma8
                Fx(:,ik1) =  exp(x(ik1))* ...
                    [Z;...
                     alpha*d0(omega(:,8))*(L*DIP);...
                     -alpha*d0(omega(:,8))*(L*DIP)];
              case 9 % sigma9
                Fx(:,ik1) =  exp(x(ik1))* ...
                    [Z;...
                     alpha*d0(omega(:,9))*(L*DIP);...
                     -alpha*d0(omega(:,9))*(L*DIP)]; 
              case 10 % kappa1
                Fx(:,ik1) = exp(x(ik1))* ...
                    [-d0(omega(:,1))*DOP;...
                     Z;...
                     d0(omega(:,1))*DOP];
              case 11 % kappa2
                Fx(:,ik1) = exp(x(ik1))* ...
                    [-d0(omega(:,2))*DOP;...
                     Z;...
                     d0(omega(:,2))*DOP];
              case 12 % kappa3
                Fx(:,ik1) = exp(x(ik1))* ...
                    [-d0(omega(:,3))*DOP;...
                     Z;...
                     d0(omega(:,3))*DOP];
              case 13 % kappa4
                Fx(:,ik1) = exp(x(ik1))* ...
                    [-d0(omega(:,4))*DOP;...
                     Z;...
                     d0(omega(:,4))*DOP];
              case 14 % kappa5
                Fx(:,ik1) = exp(x(ik1))* ...
                    [-d0(omega(:,5))*DOP;...
                     Z;...
                     d0(omega(:,5))*DOP];
              case 15 % kappa6
                Fx(:,ik1) = exp(x(ik1))* ...
                    [-d0(omega(:,6))*DOP;...
                     Z;...
                     d0(omega(:,6))*DOP];
              case 16 % kappa7
                Fx(:,ik1) = exp(x(ik1))* ...
                    [-d0(omega(:,7))*DOP;...
                     Z;...
                     d0(omega(:,7))*DOP];
              case 17 % kappa8
                Fx(:,ik1) = exp(x(ik1))* ...
                    [-d0(omega(:,8))*DOP;...
                     Z;...
                     d0(omega(:,8))*DOP];
              case 18 % kappa9
                Fx(:,ik1) = exp(x(ik1))* ...
                    [-d0(omega(:,9))*DOP;...
                     Z;...
                     d0(omega(:,9))*DOP];                
              case 19 % km1
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     -(v*d0(omega(:,1))*DOP)./Pi;...
                     (v*d0(omega(:,1))*DOP)./Pi];
              case 20 % km2
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     -(v*d0(omega(:,2))*DOP)./Pi;...
                     (v*d0(omega(:,2))*DOP)./Pi];
              case 21 % km3
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     -(v*d0(omega(:,3))*DOP)./Pi;...
                     (v*d0(omega(:,3))*DOP)./Pi];
              case 22 % km4
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     -(v*d0(omega(:,4))*DOP)./Pi;...
                     (v*d0(omega(:,4))*DOP)./Pi];
              case 23 % km5
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     -(v*d0(omega(:,5))*DOP)./Pi;...
                     (v*d0(omega(:,5))*DOP)./Pi];
              case 24 % km6
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     -(v*d0(omega(:,6))*DOP)./Pi;...
                     (v*d0(omega(:,6))*DOP)./Pi];
               case 25 % km7
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     -(v*d0(omega(:,7))*DOP)./Pi;...
                     (v*d0(omega(:,7))*DOP)./Pi];
               case 26 % km8
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     -(v*d0(omega(:,8))*DOP)./Pi;...
                     (v*d0(omega(:,8))*DOP)./Pi];
              case 27 % km9
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     -(v*d0(omega(:,9))*DOP)./Pi;...
                     (v*d0(omega(:,9))*DOP)./Pi]; 
              case 28 % alpha
                Fx(:,ik1) = exp(x(ik1))* ...
                    [L*DIP;...
                     -(1-sigma1)*L*DIP;...
                     -sigma1*L*DIP];
              case 29 %beta
                dLambdadbeta = 0*Lambda;
                dLambdadbeta(:,:,1) = log(npp).*LAM(:,:,1);
                dLambdadbeta(:,:,2) = log(npp).*LAM(:,:,2);
                iz = find(isinf(dLambdadbeta(:)));
                dLambdadbeta(iz) = 0;
                inan = find(isnan(dLambdadbeta(:)));
                dLambdadbeta(inan) = 0;
                dLdbeta = d0(dLambdadbeta(iwet)); 
                Fx(:,ik1) = exp(x(ik1))*[ alpha*dLdbeta*DIP;...
                                    -(1-sigma1)*alpha*dLdbeta*DIP;...
                                    -sigma1*alpha*dLdbeta*DIP];
                % will need L and dLdbeta for gradients of other
                % biogeochemical cycles
                parm.L = L; 
                parm.dLdbeta = dLdbeta;
                
              case 30 % bp                
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     dPFDdb*POP;...
                     Z];
            end
        end
        % Px is the derivative of the solution wrt to the parameters
        Px = mfactor(FFp,-Fx);
    end
    
    
    
    if (nargout > 2)
        %
        %
        % Compute the 2nd derivative of the solution wrt the parameters
        %
        %
        % initialize the 2nd derivative of the Jacobian wrt the parameters
        for rr = 1:3 % 3 is the number of species DIP,POP,DOP
            for cc = 1:3 % 3 is the number of species DIP,POP,DOP
                for i1 = 1:30 % 30 is the number of parameters, i.e. sigma1-9,kappa1-9,km1-9, alpha, beta, b
                    for i2 = i1:30 % 30 is the number of parameters
                        d2J{rr,cc,i1,i2} = sparse(nwet,nwet);
                    end
                end
            end
        end
        DIPx = Px(1:nwet,:);
        POPx = Px(nwet+1:2*nwet,:);
        DOPx = Px(2*nwet+1:end,:);
        % compute only the upper triangular part of the matrix
        Z = zeros(nwet,length(ip));
        for i1 = 1:length(ip)
            switch ip(i1)

              case 1 % sigma1
                FpxPx(:,:,i1) =  exp(x(i1))*[  Z;...
                                    alpha*d0(omega(:, 1))*(L*DIPx);...
                                    -alpha*d0(omega(:, 1))*(L*DIPx)];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 1 % sigma1,sigma1
                      case 2 % sigma1,sigma2
                      case 3 % sigma1,sigma3
                      case 4 % sigma1,sigma4
                      case 5 % sigma1,sigma5
                      case 6 % sigma1,sigma6
                      case 7 % sigma1,sigma7
                      case 8 % sigma1,sigma8
                      case 9 % sigma1,sigma9
                      case 10 % sigma1,kappa1
                      case 11 % sigma1,kappa2
                      case 12 % sigma1,kappa3
                      case 13 % sigma1,kappa4
                      case 14 % sigma1,kappa5
                      case 15 % sigma1,kappa6
                      case 16 % sigma1,kappa7
                      case 17 % sigma1,kappa8
                      case 18 % sigma1,kappa9
                      case 19 % sigma1,km1
                      case 20 % sigma1,km2
                      case 21 % sigma1,km3
                      case 22 % sigma1,km4
                      case 23 % sigma1,km5
                      case 24 % sigma1,km6
                      case 25 % sigma1,km7
                      case 26 % sigma1,km8
                      case 27 % sigma1,km9
                      case 28 % sigma1,alpha
                        d2J{2,1,1,28} =  L;
                        d2J{3,1,1,28} = -L;
                      case 29 % sigma1,beta
                        d2J{2,1,1,29} =  alpha*dLdbeta;
                        d2J{3,1,1,29} = -alpha*dLdbeta;
                      case 30 %sigma1,bp
                    end
                end
                
            case 2 % sigma2
                FpxPx(:,:,i1) =  exp(x(i1))*[  Z;...
                                    alpha*d0(omega(:, 2))*(L*DIPx);...
                                   -alpha*d0(omega(:, 2))*(L*DIPx)];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 2 % sigma2,sigma2
                      case 3 % sigma2,sigma3
                      case 4 % sigma2,sigma4
                      case 5 % sigma2,sigma5
                      case 6 % sigma2,sigma6
                      case 7 % sigma2,sigma7
                      case 8 % sigma2,sigma8
                      case 9 % sigma2,sigma9
                      case 10 % sigma2,kappa1
                      case 11 % sigma2,kappa2
                      case 12 % sigma2,kappa3
                      case 13 % sigma2,kappa4
                      case 14 % sigma2,kappa5
                      case 15 % sigma2,kappa6
                      case 16 % sigma2,kappa7
                      case 17 % sigma2,kappa8
                      case 18 % sigma2,kappa9
                      case 19 % sigma2,km1
                      case 20 % sigma2,km2
                      case 21 % sigma2,km3
                      case 22 % sigma2,km4
                      case 23 % sigma2,km5
                      case 24 % sigma2,km6
                      case 25 % sigma2,km7
                      case 26 % sigma2,km8
                      case 27 % sigma2,km9
                      case 28 % sigma1,alpha
                        d2J{2,1,2,28} =  L;
                        d2J{3,1,2,28} = -L;
                      case 29 % sigma1,beta
                        d2J{2,1,2,29} =  alpha*dLdbeta;
                        d2J{3,1,2,29} = -alpha*dLdbeta;
                      case 30 %sigma1,bp
                    end
                end
                
              case 3 % sigma3
                FpxPx(:,:,i1) =  exp(x(i1))*[  Z;...
                                    alpha*d0(omega(:, 3))*(L*DIPx);...
                                    -alpha*d0(omega(:, 3))*(L*DIPx)];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 3 % sigma3,sigma3
                      case 4 % sigma3,sigma4
                      case 5 % sigma3,sigma5
                      case 6 % sigma3,sigma6
                      case 7 % sigma3,sigma7
                      case 8 % sigma3,sigma8
                      case 9 % sigma3,sigma9
                      case 10 % sigma3,kappa1
                      case 11 % sigma3,kappa2
                      case 12 % sigma3,kappa3
                      case 13 % sigma3,kappa4
                      case 14 % sigma3,kappa5
                      case 15 % sigma3,kappa6
                      case 16 % sigma3,kappa7
                      case 17 % sigma3,kappa8
                      case 18 % sigma3,kappa9
                      case 19 % sigma3,km1
                      case 20 % sigma3,km2
                      case 21 % sigma3,km3
                      case 22 % sigma3,km4
                      case 23 % sigma3,km5
                      case 24 % sigma3,km6
                      case 25 % sigma3,km7
                      case 26 % sigma3,km8
                      case 27 % sigma3,km9
                      case 28 % sigma1,alpha
                        d2J{2,1,3,28} =  L;
                        d2J{3,1,3,28} = -L;
                      case 29 % sigma1,beta
                        d2J{2,1,3,29} =  alpha*dLdbeta;
                        d2J{3,1,3,29} = -alpha*dLdbeta;
                      case 30 %sigma1,bp
                    end
                end                
                
              case 4 % sigma4
                FpxPx(:,:,i1) =  exp(x(i1))*[  Z;...
                                    alpha*d0(omega(:, 4))*(L*DIPx);...
                                    -alpha*d0(omega(:, 4))*(L*DIPx)];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 4 % sigma4,sigma4
                      case 5 % sigma4,sigma5
                      case 6 % sigma4,sigma6
                      case 7 % sigma4,sigma7
                      case 8 % sigma4,sigma8
                      case 9 % sigma4,sigma9
                      case 10 % sigma4,kappa1
                      case 11 % sigma4,kappa2
                      case 12 % sigma4,kappa3
                      case 13 % sigma4,kappa4
                      case 14 % sigma4,kappa5
                      case 15 % sigma4,kappa6
                      case 16 % sigma4,kappa7
                      case 17 % sigma4,kappa8
                      case 18 % sigma4,kappa9
                      case 19 % sigma4,km1
                      case 20 % sigma4,km2
                      case 21 % sigma4,km3
                      case 22 % sigma4,km4
                      case 23 % sigma4,km5
                      case 24 % sigma4,km6
                      case 25 % sigma4,km7
                      case 26 % sigma4,km8
                      case 27 % sigma4,km9
                      case 28 % sigma1,alpha
                        d2J{2,1,4,28} =  L;
                        d2J{3,1,4,28} = -L;
                      case 29 % sigma1,beta
                        d2J{2,1,4,29} =  alpha*dLdbeta;
                        d2J{3,1,4,29} = -alpha*dLdbeta;
                      case 30 %sigma1,bp
                    end
                end
                
              case 5 % sigma5
                FpxPx(:,:,i1) =  exp(x(i1))*[  Z;...
                                    alpha*d0(omega(:, 5))*(L*DIPx);...
                                    -alpha*d0(omega(:, 5))*(L*DIPx)];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 5 % sigma5,sigma5
                      case 6 % sigma5,sigma6
                      case 7 % sigma5,sigma7
                      case 8 % sigma5,sigma8
                      case 9 % sigma5,sigma9
                      case 10 % sigma5,kappa1
                      case 11 % sigma5,kappa2
                      case 12 % sigma5,kappa3
                      case 13 % sigma5,kappa4
                      case 14 % sigma5,kappa5
                      case 15 % sigma5,kappa6
                      case 16 % sigma5,kappa7
                      case 17 % sigma5,kappa8
                      case 18 % sigma5,kappa9
                      case 19 % sigma5,km1
                      case 20 % sigma5,km2
                      case 21 % sigma5,km3
                      case 22 % sigma5,km4
                      case 23 % sigma5,km5
                      case 24 % sigma5,km6
                      case 25 % sigma5,km7
                      case 26 % sigma5,km8
                      case 27 % sigma5,km9
                      case 28 % sigma1,alpha
                        d2J{2,1,5,28} =  L;
                        d2J{3,1,5,28} = -L;
                      case 29 % sigma1,beta
                        d2J{2,1,5,29} =  alpha*dLdbeta;
                        d2J{3,1,5,29} = -alpha*dLdbeta;
                      case 30 %sigma1,bp
                    end
                end
                
              case 6 % sigma6
                FpxPx(:,:,i1) =  exp(x(i1))*[  Z;...
                                    alpha*d0(omega(:, 6))*(L*DIPx);...
                                    -alpha*d0(omega(:, 6))*(L*DIPx)];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 6 % sigma6,sigma6
                      case 7 % sigma6,sigma7
                      case 8 % sigma6,sigma8
                      case 9 % sigma6,sigma9
                      case 10 % sigma6,kappa1
                      case 11 % sigma6,kappa2
                      case 12 % sigma6,kappa3
                      case 13 % sigma6,kappa4
                      case 14 % sigma6,kappa5
                      case 15 % sigma6,kappa6
                      case 16 % sigma6,kappa7
                      case 17 % sigma6,kappa8
                      case 18 % sigma6,kappa9
                      case 19 % sigma6,km1
                      case 20 % sigma6,km2
                      case 21 % sigma6,km3
                      case 22 % sigma6,km4
                      case 23 % sigma6,km5
                      case 24 % sigma6,km6
                      case 25 % sigma6,km7
                      case 26 % sigma6,km8
                      case 27 % sigma6,km9
                      case 28 % sigma1,alpha
                        d2J{2,1,6,28} =  L;
                        d2J{3,1,6,28} = -L;
                      case 29 % sigma1,beta
                        d2J{2,1,6,29} =  alpha*dLdbeta;
                        d2J{3,1,6,29} = -alpha*dLdbeta;
                      case 30 %sigma1,bp
                    end
                end
                
              case 7 % sigma7
                FpxPx(:,:,i1) =  exp(x(i1))*[  Z;...
                                    alpha*d0(omega(:, 7))*(L*DIPx);...
                                    -alpha*d0(omega(:, 7))*(L*DIPx)];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 7 % sigma7,sigma7
                      case 8 % sigma7,sigma8
                      case 9 % sigma7,sigma9
                      case 10 % sigma7,kappa1
                      case 11 % sigma7,kappa2
                      case 12 % sigma7,kappa3
                      case 13 % sigma7,kappa4
                      case 14 % sigma7,kappa5
                      case 15 % sigma7,kappa6
                      case 16 % sigma7,kappa7
                      case 17 % sigma7,kappa8
                      case 18 % sigma7,kappa9
                      case 19 % sigma7,km1
                      case 20 % sigma7,km2
                      case 21 % sigma7,km3
                      case 22 % sigma7,km4
                      case 23 % sigma7,km5
                      case 24 % sigma7,km6
                      case 25 % sigma7,km7
                      case 26 % sigma7,km8
                      case 27 % sigma7,km9
                      case 28 % sigma1,alpha
                        d2J{2,1,7,28} =  L;
                        d2J{3,1,7,28} = -L;
                      case 29 % sigma1,beta
                        d2J{2,1,7,29} =  alpha*dLdbeta;
                        d2J{3,1,7,29} = -alpha*dLdbeta;
                      case 30 %sigma1,bp
                    end
                end
                
              case 8 % sigma8
                FpxPx(:,:,i1) =  exp(x(i1))*[  Z;...
                                    alpha*d0(omega(:, 8))*(L*DIPx);...
                                    -alpha*d0(omega(:, 8))*(L*DIPx)];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 8 % sigma8,sigma8
                      case 9 % sigma8,sigma9
                      case 10 % sigma8,kappa1
                      case 11 % sigma8,kappa2
                      case 12 % sigma8,kappa3
                      case 13 % sigma8,kappa4
                      case 14 % sigma8,kappa5
                      case 15 % sigma8,kappa6
                      case 16 % sigma8,kappa7
                      case 17 % sigma8,kappa8
                      case 18 % sigma8,kappa9
                      case 19 % sigma8,km1
                      case 20 % sigma8,km2
                      case 21 % sigma8,km3
                      case 22 % sigma8,km4
                      case 23 % sigma8,km5
                      case 24 % sigma8,km6
                      case 25 % sigma8,km7
                      case 26 % sigma8,km8
                      case 27 % sigma8,km9
                      case 28 % sigma1,alpha
                        d2J{2,1,8,28} =  L;
                        d2J{3,1,8,28} = -L;
                      case 29 % sigma1,beta
                        d2J{2,1,8,29} =  alpha*dLdbeta;
                        d2J{3,1,8,29} = -alpha*dLdbeta;
                      case 30 %sigma1,bp
                    end
                end
                
              case 9 % sigma9
                FpxPx(:,:,i1) =  exp(x(i1))*[  Z;...
                                    alpha*d0(omega(:, 9))*(L*DIPx);...
                                    -alpha*d0(omega(:, 9))*(L*DIPx)];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 9 % sigma9,sigma9
                      case 10 % sigma9,kappa1
                      case 11 % sigma9,kappa2
                      case 12 % sigma9,kappa3
                      case 13 % sigma9,kappa4
                      case 14 % sigma9,kappa5
                      case 15 % sigma9,kappa6
                      case 16 % sigma9,kappa7
                      case 17 % sigma9,kappa8
                      case 18 % sigma9,kappa9
                      case 19 % sigma9,km1
                      case 20 % sigma9,km2
                      case 21 % sigma9,km3
                      case 22 % sigma9,km4
                      case 23 % sigma9,km5
                      case 24 % sigma9,km6
                      case 25 % sigma9,km7
                      case 26 % sigma9,km8
                      case 27 % sigma9,km9
                      case 28 % sigma1,alpha
                        d2J{2,1,9,28} =  L;
                        d2J{3,1,9,28} = -L;
                      case 29 % sigma1,beta
                        d2J{2,1,9,29} =  alpha*dLdbeta;
                        d2J{3,1,9,29} = -alpha*dLdbeta;
                      case 30 %sigma1,bp
                    end
                end
                
              case 10 % kappa1
                FpxPx(:,:,i1) = exp(x(i1))*[-d0(omega(:, 1))*DOPx;...
                                    Z;...
                                    d0(omega(:, 1))*DOPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 10 % kappa1,kappa1
                      case 11 % kappa1,kappa2
                      case 12 % kappa1,kappa3
                      case 13 % kappa1,kappa4
                      case 14 % kappa1,kappa5
                      case 15 % kappa1,kappa6
                      case 16 % kappa1,kappa7
                      case 17 % kappa1,kappa8
                      case 18 % kappa1,kappa9
                      case 19 % kappa1,km1
                      case 20 % kappa1,km2
                      case 21 % kappa1,km3
                      case 22 % kappa1,km4
                      case 23 % kappa1,km5
                      case 24 % kappa1,km6
                      case 25 % kappa1,km7
                      case 26 % kappa1,km8
                      case 27 % kappa1,km9
                      case 28 % kappa1,alpha
                      case 29 % kappa1,beta
                      case 30 % kappa1,bp
                    end
                end
                
              case 11 % kappa2
                FpxPx(:,:,i1) = exp(x(i1))*[-d0(omega(:, 2))*DOPx;...
                                    Z;...
                                    d0(omega(:, 2))*DOPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 11 % kappa2,kappa2
                      case 12 % kappa2,kappa3
                      case 13 % kappa2,kappa4
                      case 14 % kappa2,kappa5
                      case 15 % kappa2,kappa6
                      case 16 % kappa2,kappa7
                      case 17 % kappa2,kappa8
                      case 18 % kappa2,kappa9
                      case 19 % kappa2,km1
                      case 20 % kappa2,km2
                      case 21 % kappa2,km3
                      case 22 % kappa2,km4
                      case 23 % kappa2,km5
                      case 24 % kappa2,km6
                      case 25 % kappa2,km7
                      case 26 % kappa2,km8
                      case 27 % kappa2,km9
                      case 28 % kappa1,alpha
                      case 29 % kappa1,beta
                      case 30 % kappa1,bp
                    end
                end
                
              case 12 % kappa3
                FpxPx(:,:,i1) = exp(x(i1))*[-d0(omega(:, 3))*DOPx;...
                                    Z;...
                                    d0(omega(:, 3))*DOPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 12 % kappa3,kappa3
                      case 13 % kappa3,kappa4
                      case 14 % kappa3,kappa5
                      case 15 % kappa3,kappa6
                      case 16 % kappa3,kappa7
                      case 17 % kappa3,kappa8
                      case 18 % kappa3,kappa9
                      case 19 % kappa3,km1
                      case 20 % kappa3,km2
                      case 21 % kappa3,km3
                      case 22 % kappa3,km4
                      case 23 % kappa3,km5
                      case 24 % kappa3,km6
                      case 25 % kappa3,km7
                      case 26 % kappa3,km8
                      case 27 % kappa3,km9
                      case 28 % kappa1,alpha
                      case 29 % kappa1,beta
                      case 30 % kappa1,bp
                    end
                end                
                
              case 13 % kappa4
                FpxPx(:,:,i1) = exp(x(i1))*[-d0(omega(:, 4))*DOPx;...
                                    Z;...
                                    d0(omega(:, 4))*DOPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 13 % kappa4,kappa4
                      case 14 % kappa4,kappa5
                      case 15 % kappa4,kappa6
                      case 16 % kappa4,kappa7
                      case 17 % kappa4,kappa8
                      case 18 % kappa4,kappa9
                      case 19 % kappa4,km1
                      case 20 % kappa4,km2
                      case 21 % kappa4,km3
                      case 22 % kappa4,km4
                      case 23 % kappa4,km5
                      case 24 % kappa4,km6
                      case 25 % kappa4,km7
                      case 26 % kappa4,km8
                      case 27 % kappa4,km9
                      case 28 % kappa1,alpha
                      case 29 % kappa1,beta
                      case 30 % kappa1,bp
                    end
                end 
                
              case 14 % kappa5
                FpxPx(:,:,i1) = exp(x(i1))*[-d0(omega(:, 5))*DOPx;...
                                    Z;...
                                    d0(omega(:, 5))*DOPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 14 % kappa5,kappa5
                      case 15 % kappa5,kappa6
                      case 16 % kappa5,kappa7
                      case 17 % kappa5,kappa8
                      case 18 % kappa5,kappa9
                      case 19 % kappa5,km1
                      case 20 % kappa5,km2
                      case 21 % kappa5,km3
                      case 22 % kappa5,km4
                      case 23 % kappa5,km5
                      case 24 % kappa5,km6
                      case 25 % kappa5,km7
                      case 26 % kappa5,km8
                      case 27 % kappa5,km9
                      case 28 % kappa1,alpha
                      case 29 % kappa1,beta
                      case 30 % kappa1,bp
                    end
                end  
                
              case 15 % kappa6
                FpxPx(:,:,i1) = exp(x(i1))*[-d0(omega(:, 6))*DOPx;...
                                    Z;...
                                    d0(omega(:, 6))*DOPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 15 % kappa6,kappa6
                      case 16 % kappa6,kappa7
                      case 17 % kappa6,kappa8
                      case 18 % kappa6,kappa9
                      case 19 % kappa6,km1
                      case 20 % kappa6,km2
                      case 21 % kappa6,km3
                      case 22 % kappa6,km4
                      case 23 % kappa6,km5
                      case 24 % kappa6,km6
                      case 25 % kappa6,km7
                      case 26 % kappa6,km8
                      case 27 % kappa6,km9
                      case 28 % kappa1,alpha
                      case 29 % kappa1,beta
                      case 30 % kappa1,bp
                    end
                end  
                
              case 16 % kappa7
                FpxPx(:,:,i1) = exp(x(i1))*[-d0(omega(:, 7))*DOPx;...
                                    Z;...
                                    d0(omega(:, 7))*DOPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 16 % kappa7,kappa7
                      case 17 % kappa7,kappa8
                      case 18 % kappa7,kappa9
                      case 19 % kappa7,km1
                      case 20 % kappa7,km2
                      case 21 % kappa7,km3
                      case 22 % kappa7,km4
                      case 23 % kappa7,km5
                      case 24 % kappa7,km6
                      case 25 % kappa7,km7
                      case 26 % kappa7,km8
                      case 27 % kappa7,km9
                      case 28 % kappa1,alpha
                      case 29 % kappa1,beta
                      case 30 % kappa1,bp
                    end
                end   
                
              case 17 % kappa8
                FpxPx(:,:,i1) = exp(x(i1))*[-d0(omega(:, 8))*DOPx;...
                                    Z;...
                                    d0(omega(:, 8))*DOPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 17 % kappa8,kappa8
                      case 18 % kappa8,kappa9
                      case 19 % kappa8,km1
                      case 20 % kappa8,km2
                      case 21 % kappa8,km3
                      case 22 % kappa8,km4
                      case 23 % kappa8,km5
                      case 24 % kappa8,km6
                      case 25 % kappa8,km7
                      case 26 % kappa8,km8
                      case 27 % kappa8,km9
                      case 28 % kappa1,alpha
                      case 29 % kappa1,beta
                      case 30 % kappa1,bp
                    end
                end  
                
              case 18 % kappa9
                FpxPx(:,:,i1) = exp(x(i1))*[-d0(omega(:, 9))*DOPx;...
                                    Z;...
                                    d0(omega(:, 9))*DOPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 18 % kappa9,kappa9
                      case 19 % kappa9,km1
                      case 20 % kappa9,km2
                      case 21 % kappa9,km3
                      case 22 % kappa9,km4
                      case 23 % kappa9,km5
                      case 24 % kappa9,km6
                      case 25 % kappa9,km7
                      case 26 % kappa9,km8
                      case 27 % kappa9,km9
                      case 28 % kappa1,alpha
                      case 29 % kappa1,beta
                      case 30 % kappa1,bp
                    end
                end               
                
              case 19 % km1
                FpxPx(:,:,i1) = exp(x(i1))*[Z;...
                                    -v*d0(omega(:, 1))*DOPx./Pi;...
                                    v*d0(omega(:, 1))*DOPx./Pi];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 19 % km1,km1
                      case 20 % km1,km2
                      case 21 % km1,km3
                      case 22 % km1,km4
                      case 23 % km1,km5
                      case 24 % km1,km6
                      case 25 % km1,km7
                      case 26 % km1,km8
                      case 27 % km1,km9
                      case 28 % km1,alpha
                      case 29 % km1,beta
                      case 30 % km1,bp
                    end
                end
                
              case 20 % km2
                FpxPx(:,:,i1) = exp(x(i1))*[Z;...
                                    -v*d0(omega(:, 2))*DOPx./Pi;...
                                    v*d0(omega(:, 2))*DOPx./Pi];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 20 % km2,km2
                      case 21 % km2,km3
                      case 22 % km2,km4
                      case 23 % km2,km5
                      case 24 % km2,km6
                      case 25 % km2,km7
                      case 26 % km2,km8
                      case 27 % km2,km9
                      case 28 % km1,alpha
                      case 29 % km1,beta
                      case 30 % km1,bp
                    end
                end 
                
              case 21 % km3
                FpxPx(:,:,i1) = exp(x(i1))*[Z;...
                                    -v*d0(omega(:, 3))*DOPx./Pi;...
                                    v*d0(omega(:, 3))*DOPx./Pi];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 21 % km3,km3
                      case 22 % km3,km4
                      case 23 % km3,km5
                      case 24 % km3,km6
                      case 25 % km3,km7
                      case 26 % km3,km8
                      case 27 % km3,km9
                      case 28 % km1,alpha
                      case 29 % km1,beta
                      case 30 % km1,bp
                    end
                end   
                
              case 22 % km4
                FpxPx(:,:,i1) = exp(x(i1))*[Z;...
                                    -v*d0(omega(:, 4))*DOPx./Pi;...
                                    v*d0(omega(:, 4))*DOPx./Pi];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 22 % km4,km4
                      case 23 % km4,km5
                      case 24 % km4,km6
                      case 25 % km4,km7
                      case 26 % km4,km8
                      case 27 % km4,km9
                      case 28 % km1,alpha
                      case 29 % km1,beta
                      case 30 % km1,bp
                    end
                end   
                
              case 23 % km5
                FpxPx(:,:,i1) = exp(x(i1))*[Z;...
                                    -v*d0(omega(:, 5))*DOPx./Pi;...
                                    v*d0(omega(:, 5))*DOPx./Pi];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 23 % km5,km5
                      case 24 % km5,km6
                      case 25 % km5,km7
                      case 26 % km5,km8
                      case 27 % km5,km9
                      case 28 % km1,alpha
                      case 29 % km1,beta
                      case 30 % km1,bp
                    end
                end  
                
              case 24 % km6
                FpxPx(:,:,i1) = exp(x(i1))*[Z;...
                                    -v*d0(omega(:, 6))*DOPx./Pi;...
                                    v*d0(omega(:, 6))*DOPx./Pi];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 24 % km6,km6
                      case 25 % km6,km7
                      case 26 % km6,km8
                      case 27 % km6,km9
                      case 28 % km1,alpha
                      case 29 % km1,beta
                      case 30 % km1,bp
                    end
                end 
                
              case 25 % km7
                FpxPx(:,:,i1) = exp(x(i1))*[Z;...
                                    -v*d0(omega(:, 7))*DOPx./Pi;...
                                    v*d0(omega(:, 7))*DOPx./Pi];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 25 % km7,km7
                      case 26 % km7,km8
                      case 27 % km7,km9
                      case 28 % km1,alpha
                      case 29 % km1,beta
                      case 30 % km1,bp
                    end
                end 
                
              case 26 % km8
                FpxPx(:,:,i1) = exp(x(i1))*[Z;...
                                    -v*d0(omega(:, 8))*DOPx./Pi;...
                                    v*d0(omega(:, 8))*DOPx./Pi];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 26 % km8,km8
                      case 27 % km8,km9
                      case 28 % km1,alpha
                      case 29 % km1,beta
                      case 30 % km1,bp
                    end
                end  
                
              case 27 % km9
                FpxPx(:,:,i1) = exp(x(i1))*[Z;...
                                    -v*d0(omega(:, 9))*DOPx./Pi;...
                                    v*d0(omega(:, 9))*DOPx./Pi];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 27 % km9,km9
                      case 28 % km1,alpha
                      case 29 % km1,beta
                      case 30 % km1,bp
                    end
                end
                          
               case 28 % alpha
                 FpxPx(:,:,i1) = exp(x(i1))*[  L*DIPx;...
                                    -(1-sigma1)*L*DIPx;...
                                    -sigma1*L*DIPx];   
                 for i2 = i1:length(ip)
                    switch ip(i2)
                      case 28 % alpha,alpha
                      case 29 % alpha,beta
                        d2J{1,1,28,29} = dLdbeta;
                        d2J{2,1,28,29} = -(1-sigma1)*dLdbeta;
                        d2J{3,1,28,29} =     -sigma1*dLdbeta;
                      case 30 % alpha,bp
                    end
                 end
                 
            case 29 % beta
                FpxPx(:,:,i1) = exp(x(i1))*[  alpha*dLdbeta*DIPx;...
                                    -(1-sigma1)*alpha*dLdbeta*DIPx;...
                                    -sigma1*alpha*dLdbeta*DIPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 29 % beta,beta
                        d2Lambdadbetadbeta = 0*Lambda;
                        d2Lambdadbetadbeta(:,:,1) = log(npp).*log(npp).*LAM(:,:,1);
                        d2Lambdadbetadbeta(:,:,2) = log(npp).*log(npp).*LAM(:,:,2);
                        iz = find(isinf(d2Lambdadbetadbeta(:)));
                        d2Lambdadbetadbeta(iz) = 0;
                        inan = find(isnan(d2Lambdadbetadbeta(:)));
                        d2Lambdadbetadbeta(inan) = 0;
                        d2Ldbetadbeta = d0(d2Lambdadbetadbeta(iwet)); 
                        d2J{1,1,29,29} = alpha*d2Ldbetadbeta;
                        d2J{2,1,29,29} = -(1-sigma1)*alpha*d2Ldbetadbeta;
                        d2J{3,1,29,29} =     -sigma1*alpha*d2Ldbetadbeta;
                        parm.d2Ldbetadbeta = d2Ldbetadbeta;
                      case 30 % beta,bp
                    end
                end
                
              case 30 % bp
                FpxPx(:,:,i1) = exp(x(i1))*[  Z;...
                                    dPFDdb*POPx;...
                                    Z];
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 30 % bp,bp
                        d2J{2,2,30,30} = d2PFDdb2;
                    end
                end               
            end
        end
        symcol = @(i,j) (i>=j).*i+(i<j).*j;
        symrow = @(i,j) (i>=j).*j+(i<j).*i;
        delta = @(i,j) (i==j)*1 + (i~=j)*0.0;
        k = 0;
        for i1 = 1:length(ip)
            for i2 = i1:length(ip) 
                k = k+1;
                rhs(:,k) = -(...
                    exp(x(i1)+x(i2))*cell2mat(d2J(:,:,symrow(ip(i1),ip(i2)),symcol(ip(i1),ip(i2))))*P+...
                    FpxPx(:,i1,i2)+FpxPx(:,i2,i1)+delta(i1,i2)*Fx(:,i1));

            end
        end
        % Pxx is the second derivative of the solution wrt to the parameters
        Pxx = mfactor(FFp,rhs);
    end
end

