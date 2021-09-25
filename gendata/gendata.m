% This script generates the input files for a trial MITgcm simulation using
% the 'IcePlume' package



% Accuracy of binary files
acc = 'real*8';

% Number of time levels for time varying forcing
nt = 1;


%create grid and bathymetry

deltaX=150;
deltaY=150;
nx=round(delxold(1)/deltaX* nxold);
ny=round(delyold(1)/deltaY* nyold);



[X,Y,Z] = meshgrid(cumsum(delyold),cumsum(delxold), linspace(0,deltaZ*nz, nz) );

[Xq,Yq,Zq] = meshgrid( linspace(0,deltaY*ny, ny),linspace(0,deltaX*nx, nx),linspace(0,deltaZ*nz, nz));



%%%First create the bathymetry from the bedmachine data, separately




%do some very adhoc post-processing to the bathymetry, mostly to smooth it
%around the periodic edges and to remove edge regions outside of the
%periodic boundary

generate_Tor_bath_strat; %this will generate the raw bathymetry and OMG stratification inputs


bathy_new = newbathy(11:end,:);


%bathy_new(end-16:end-10,240:end)=zeros;%%%should smooth this out more
bathy_new(1:27,1:43)=zeros;%%%should smooth this out more
%bathy_new(200:end,1:15)=zeros;

bathy_new(:,end-10:end)=zeros;
bathy_new(:,1:10)=zeros;
bathy_new(216,86)=bathy_new(217,86);



bath=bathy_new;

 bath_new=bath;
 bath_new(bath<-200) = -200;
 bathsmooth=smooth2a(bath_new,3,3);
 
 for ii=1:size(bath,1)
     for jj=1:size(bath,2)
         if bath(ii,jj)<-200
 bathsmooth(ii,jj)=bath(ii,jj);
         end
      
     end
 end

  
 for ii=1:size(bath,1)
     for jj=1:size(bath,2)
         
          if bathsmooth(ii,jj)>-20
 bathsmooth(ii,jj)=0;
         end
     end
 end         
 
 bathsmooth(272:end,:)=bathy_new(272:end,:);
 

bathsmooth(:,1)= bathsmooth(:,2);
bathsmooth(1,:) = bathsmooth(2,:);

fid=fopen('bathymetry.bin','w','b'); fwrite(fid,bathsmooth,acc);fclose(fid);




%print out the delx and dely bins

% Dimensions of grid
nx=size(bathsmooth,1);%360+220;  %10 cores each 58 x
ny=size(bathsmooth,2);%560+98;   %14 cores each 47


% x scale 2232000
delx = zeros(1,nx);
delx(:) = deltaX;

%delx(1,end-round(nx-50000/deltaX):end) = linspace(deltaX, 8*deltaX, round(nx-50000/deltaX)+1);
lx=sum(delx);

fid=fopen('delx.bin','w','b'); fwrite(fid,delx,acc);fclose(fid);


% y scale
dely = zeros(1,ny);
dely(:) = deltaY;


ly=sum(dely);

fid=fopen('dely.bin','w','b'); fwrite(fid,dely,acc);fclose(fid);






%% Temperature and salinity profiles
% These are used to write initial conditions, boundary conditions etc
% using csv files from OMG data, which is used to specify open-ocean
% nudging as well as initial data



saltt_=[saltt; repmat(saltt(end),-deepest-size(Z_,1),1)];
tempp_=[tempp; repmat(tempp(end),-deepest-size(Z_,1),1)];

saltprof=zeros(1,Nr);
tempprof=zeros(1,Nr);
for i=1:95
    saltprof(i) = nanmean(saltt_(  round(deltaZ*(i-1)+1): round(deltaZ*(i))));
    tempprof(i) = nanmean(tempp_(  round(deltaZ*(i-1)+1): round(deltaZ*(i))));
end
 for i=96:100
   saltprof(i)=saltprof(95);
   tempprof(i)=tempprof(95);
 end
   

saltini = zeros(nx,ny,nz);
tempini = zeros(nx,ny,nz);

saltini= permute(repmat(saltprof',1,nx,ny),[2 3 1]);
tempini= permute(repmat(tempprof',1,nx,ny),[2 3 1]);

figure, subplot(1,2,1)
plot(tempprof) , hold on
plot(tempp_(1:deltaZ:end))
subplot(1,2,2)
plot(saltprof), hold on
plot(saltt_(1:deltaZ:end))

fid=fopen('saltini.bin','w','b'); fwrite(fid,saltini,acc);fclose(fid);
fid=fopen('tempini.bin','w','b'); fwrite(fid,tempini,acc);fclose(fid);

%% Subglacial runoff
% the granular data is based on Mankoff 2017 dataset and takes an average
% discharge over the summer months

% Define the velocity (m/s) of subglacial runoff as it enters the fjord.
% 1 m/s seems a reasonable value (results not sensitive to this value).
wsg = 1;

% Templates
runoffVel   = zeros(nx,ny);
runoffRad   = zeros(nx,ny);
plumeMask   = zeros(nx,ny);

%%% Define plume-type mask %%%
% 1 = ice but no plume (melting only)
% 2 = sheet plume (Jenkins 2011)
% 3 = half-conical plume (Cowton et al 2015)
% 4 = both sheet plume and half-conical plume (NOT YET IMPLEMENTED)
% 5 = detaching conical plume (Goldberg)

% POSITIVE values indicate ice front is orientated north-south
% NEGATIVE values indicate ice front is orientated east-west

% Define melting along the glacier front (located at fjord head)

%plumeMask(icefront,2:(end-1)) = 2;
%plumeMask(icefront,round(ny/2) -fjwidth/deltaY/2+npy :  round(ny/2)+fjwidth/deltaY/2-npy  ) = 0;


% The plume will be located in the fjord centre at the glacier
% front
%% I have not fully automated this process yet, which is very annoying!!! and will be fixed in a future version


%locations of melt plumes (all along glacial face)
plumeMask(560:569,231)=2;
plumeMask(569,229:231)=2;
plumeMask(569:572,229)=2;
plumeMask(572,229:231)=2;
plumeMask(572:583,231)=2;
plumeMask(583,229:231)=2;
plumeMask(583:603,229)=2;

%location of discharge plume)
plumeMask(572,229)=3;


plumeMask(619:627,128)=2;
plumeMask(627,128:130)=2;
plumeMask(627:633,130)=2;
plumeMask(633,130:132)=2;
plumeMask(633:637,132)=2;
plumeMask(637,132:134)=2;
plumeMask(637:639,134)=2;
plumeMask(639,134:136)=2;
plumeMask(639:641,136)=2;
plumeMask(641,136:140)=2;
plumeMask(641:643,140)=2;
plumeMask(643,140:146)=2;
plumeMask(643:645,146)=2;
plumeMask(645,146:150)=2;
plumeMask(645:647,150)=2;
plumeMask(647,150:151)=2;
plumeMask(645:647,151)=2;
plumeMask(645,151:153)=2;
plumeMask(643:645,153)=2;
plumeMask(643,153:162)=2;
plumeMask(643:645,162)=2;
plumeMask(645,162:163)=2;
plumeMask(643:645,163)=2;
plumeMask(643,163:168)=2;
plumeMask(643:645,168)=2;
plumeMask(645,168:176)=2;


plumeMask(643,156)=3;



% Specify runoff (m^3/s)
runoff = 220;

% Define runoff velocity in each location (as specified above)

runoffVel(572,229)= wsg;

% Calculate channel radius to give runoff at velocity
%runoffRad(icefront,round(ny/2)) = sqrt(2*runoff/(pi*wsg));

runoffRad(572,229) = sqrt(2*runoff/(pi*wsg));



%An additional runoff location for the second glacier
runoff = 130;


runoffVel(643,156)= wsg;

runoffRad(643,156) = sqrt(2*runoff/(pi*wsg));



% Write files.
fid=fopen('runoffVel.bin','w','b'); fwrite(fid,runoffVel,acc);fclose(fid);
fid=fopen('runoffRad.bin','w','b'); fwrite(fid,runoffRad,acc);fclose(fid);
fid=fopen('plumeMask.bin','w','b'); fwrite(fid,plumeMask,acc);fclose(fid);





%% Boundary conditions (T/S/U/V), based on the OMG data above.
%The U/V is based on existing shelf currents/extrapolation of SSH measures data/ECCO where available and is
%more ad-hoc. However, the main goal is to to minimize inaccuracies in shelf
%variability without creating unnecessarily high along-shelf coastal
%transport


shelfz=nz;
clear EBCu EBCs EBCt  WBCu WBCs WBCt NBCt NBCu NBCs SBCt SBCu SBCs 
EBCu = zeros(ny,shelfz);
EBCs = zeros(ny,shelfz);
EBCt = zeros(ny,shelfz);
EBCv = zeros(ny,shelfz);

WBCu = zeros(ny,shelfz);
WBCs = zeros(ny,shelfz);
WBCt = zeros(ny,shelfz);
WBCv = zeros(ny,shelfz);

NBCu = zeros(nx,shelfz);
NBCs = zeros(nx,shelfz);
NBCt = zeros(nx,shelfz);
NBCv = zeros(nx,shelfz);

SBCu = zeros(nx,shelfz);
SBCs = zeros(nx,shelfz);
SBCt = zeros(nx,shelfz);
SBCv = zeros(nx,shelfz);

for i = 1:nz %%%%%%%%%%%%%
        WBCt(:,i) = tempprof(i);
        WBCs(:,i) = saltprof(i);
        WBCv(:,i) = 0.00;
end

for i = 1:nz %%%%%%%%%%%%%
        EBCt(:,i) = tempprof(i);
        EBCs(:,i) = saltprof(i);
        EBCu(:,i) = 0;

end


for i = 1:nz %%%%%%%%%%%%%
        SBCt(:,i) = tempprof(i);
        SBCs(:,i) = saltprof(i);
          SBCv(:,i) = .00;
        
        NBCt(:,i) = tempprof(i);
        NBCs(:,i) = saltprof(i);
          NBCv(:,i) = .00;
end



% Apply barotropic velocity to balance input of runoff

%fjordMouthCrossSection = -sum(bathymetry(end,:))*deltaY;
%fjordMouthVelocity = runoff/fjordMouthCrossSection;


fid=fopen('EBCu.bin','w','b'); fwrite(fid,EBCu,acc);fclose(fid);
fid=fopen('EBCs.bin','w','b'); fwrite(fid,EBCs,acc);fclose(fid);
fid=fopen('EBCt.bin','w','b'); fwrite(fid,EBCt,acc);fclose(fid);
fid=fopen('EBCv.bin','w','b'); fwrite(fid,EBCv,acc);fclose(fid);


fid=fopen('NBCu.bin','w','b'); fwrite(fid,NBCu,acc);fclose(fid);
fid=fopen('NBCs.bin','w','b'); fwrite(fid,NBCs,acc);fclose(fid);
fid=fopen('NBCt.bin','w','b'); fwrite(fid,NBCt,acc);fclose(fid);
fid=fopen('NBCv.bin','w','b'); fwrite(fid,NBCv,acc);fclose(fid);

fid=fopen('SBCu.bin','w','b'); fwrite(fid,SBCu,acc);fclose(fid);
fid=fopen('SBCs.bin','w','b'); fwrite(fid,SBCs,acc);fclose(fid);
fid=fopen('SBCt.bin','w','b'); fwrite(fid,SBCt,acc);fclose(fid);
fid=fopen('SBCv.bin','w','b'); fwrite(fid,SBCv,acc);fclose(fid);

fid=fopen('WBCu.bin','w','b'); fwrite(fid,WBCu,acc);fclose(fid);
fid=fopen('WBCs.bin','w','b'); fwrite(fid,WBCs,acc);fclose(fid);
fid=fopen('WBCt.bin','w','b'); fwrite(fid,WBCt,acc);fclose(fid);
fid=fopen('WBCv.bin','w','b'); fwrite(fid,WBCv,acc);fclose(fid);



%ptracers, used if need tracer concentrations (not currently used)

%ptracers = zeros(nx,ny,nz);
%ptracers(:,1:10,15:end)=repmat([1:-.1:.1],[nx,1,86]);


%fid=fopen('ptracerinit.bin','w','b'); fwrite(fid,ptracers,acc);fclose(fid);

%tracernudg=ones(nx,ny,nz);

%fid=fopen('ptracernudg.bin','w','b'); fwrite(fid,tracernudg,acc);fclose(fid);

%tracermask=zeros(nx,ny,nz);
%tracermask(:,1:10,15:end)=repmat([1:-.1:.1],[nx,1,86]);
%fid=fopen('ptracermask.bin','w','b'); fwrite(fid,tracermask,acc);fclose(fid);


%%%%lagrangian tracers (not currently used)

%{
xnum=10;%
znum=10;

floats=xnum*znum;
flt=zeros(9,floats);

flt_firstline=[floats -1 0 0 -1 floats 0 0 -1]';

%{ 
npart   A unique float identifier (1,2,3,...)
tstart  start date of integration of float (in s)
  - If tstart=-1 floats are integrated right from the beginning
xpart   x position of float (in units of XC)
ypart   y position of float (in units of YC)
kpart   actual vertical level of float
kfloat  target level of float
       - should be the same as kpart at the beginning
iup     flag if the float
  - should profile   ( >  0 = return cycle (in s) to surface)
  - remain at depth  ( =  0 )
  - is a 3D float    ( = -1 ).
  - should be advected WITHOUT additional noise ( = -2 ).
    (This implies that the float is non-profiling)
  - is a mooring     ( = -3 ), i.e. the float is not advected
itop    time of float the surface (in s)
tend    end date of integration of float (in s)
  - If tend=-1 floats are integrated till the end of the integration
%}


flt(1,:) = [1:xnum*znum];% float index
flt(2,:) = 2980800*ones;% float index


flt(3,:) = reshape(repmat(linspace(5, 141,xnum),[znum,1])*150,[1 xnum*znum]);

flt(4,:) = 300*ones;
flt(5,:) = reshape(repmat(linspace(-150, -350,xnum)',[1,znum]),[1 xnum*znum]);
flt(6,:) = reshape(repmat(linspace(-150, -350,xnum)',[1,znum]),[1 xnum*znum]);
flt(7,:) = -2*ones;

flt(9,:) = -1*ones;

fltfull=[flt_firstline flt];

fid=fopen('flt_ini_pos.bin','w','b'); fwrite(fid,fltfull,acc);fclose(fid);
%}




