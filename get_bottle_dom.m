function DOP = get_bottle_dom(grid,M3d,d,dataset)
% updated by Zhou Liang 04/23/2021, fix first and last vertical level bug
% Interpolate the DOP bottle data onto the model grid
% load the observations
  D = importdata(dataset);
  yobs = D.data(:,1);
  xobs = D.data(:,2);
  zobs = D.data(:,3);
  
  dopobs = D.data(:,d);
  %eliminate negative DOM value
  dopobs(dopobs <=0 == -999);
  % eliminate the missing data points
  ikp = find(dopobs~= -999);
  xobs = (xobs(ikp)>=0).*xobs(ikp)+(xobs(ikp)<0).*(xobs(ikp)+360);
  yobs = yobs(ikp);
  zobs = zobs(ikp);
  dopobs = dopobs(ikp);
  
  %
  % bin the data in the vertical
  % 
  
  % give the model grid boxes a vertical index
  zt = squeeze(grid.ZT3d(1,1,:));
  zndx = 1:length(zt);
  
  % interp the vertical model index onto the bottle coordinates
  obs_zndx = 0*dopobs;
  fprintf('Vertical indexing...'); tic
    obs_zndx = interp1(zt,zndx,zobs,'nearest');
  toc
  
  obs_zndx(find(zobs<18.1))   = 1;  % assign surface data to level 1
  obs_zndx(find(zobs>5433.3)) = 24; % assign deep ocean data to level 24
  
  
  % give the model grid boxes a horizontal index
  YT = grid.YT3d(:,:,1);
  XT = grid.XT3d(:,:,1);
  hndx = 0*XT;
  hndx(:) = 1:length(XT(:));
  
  % interp the horizontal model index onto the bottle coordinates
  obs_hndx = 0*dopobs;
  fprintf('Horizontal indexing...'); tic
    F = scatteredInterpolant(XT(:),YT(:),hndx(:),'nearest')
    obs_hndx = F(xobs(:),yobs(:));
  toc
  %
  [ny,nx,nz] = size(M3d);
  n = ny*nx*nz;
  
  % assign to each observation a unique model grid box index
  Zndx = zeros(1,1,nz);
  Zndx(1,1,:) = zndx;
  mxhndx = max(hndx(:));
  obsndx = obs_hndx+(obs_zndx-1)*mxhndx;
  [srt_ndx,indx] = sort(obsndx);
  ikp = find(~isnan(srt_ndx));
  
  % take advantage of the fact that "sparse" adds elements with duplicate
  % entries of i and j
  DOP = M3d;
  DOP(:) = full(diag(sparse(srt_ndx(ikp),srt_ndx(ikp),dopobs(indx(ikp)),n,n)));
  N = M3d;
  N(:) = full(diag(sparse(srt_ndx(ikp),srt_ndx(ikp),0*dopobs(indx(ikp))+1,n,n)));
  
  DOP = DOP./N;
  
  