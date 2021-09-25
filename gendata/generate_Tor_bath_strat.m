
%%%%%read greenland bathymetry and stratification from bedmachine and OMG
%%%%%data


%%

%netcdf file 
filename = '/Volumes/Lacie/Greenland/matlab/BedMachineGreenland-2017-09-20.nc'; %%%read in the bedmachine file
x = ncread(filename,'x');
y = ncread(filename,'y');
mask = ncread(filename,'mask')';
bed = ncread(filename,'bed')'; %Do not forget to transpose (MATLAB is column oriented)
icethickness= ncread(filename,'thickness')';
surface= ncread(filename,'surface')';
%Display bed elevation
%imagesc(x,y,bed); axis xy equal; caxis([-500 1800]);


%mask with above sea level ignored
bedmask=bed.*(mask==0);

icethickness=double(icethickness);
icesurf=(icethickness+bed).*(mask==2);




%%%%%%zoom in on torsukataq
%mask with above sea level ignored
bedmask=bed.*(mask==0);

icethickness=double(icethickness);
icesurf=(icethickness+bed).*(mask==2);

xl=2540+158; xr=3150-(610-492-15);  % from left to right
yl=10160+(871-764-15); yr=11030-662+15; %from top to bottom
imagesc(x(xl:xr), y(yl:yr), bed(yl:yr, xl:xr));
%hold on, contour(x(xl:xr),y(yl:yr), flipud(icethickness(yl:yr, xl:xr)),[0 1],'r','linewidth',3);
colormap(haxby)
caxis([-800 0]), axis xy equal;
makepretty

figure

newbathy=flipud(bed(yl:yr,xl:xr).*(mask(yl:yr,xl:xr)==0))';
pcolor(newbathy') , shading interp, colorbar
colormap(haxby)

deepest=min(min(newbathy)); %%%deepest part of the domain





%%%%%%torsukataq stratification


%%%read in T/S profiles for corresponding OMG data
%[tempp saltt Z_]=csvimport('~/Downloads/CTD_20200825_1457.csv', 'columns', [4,5,6],'noHeader',true); 
[tempp saltt Z_]=csvimport('/Volumes/Elements/OMGdata/Data/Both/CTD_20200825_1457.csv', 'columns', [4,5,6],'noHeader',true); 
figure
plot(tempp, [-1:-1: -737]), hold on
plot(smooth(tempp(1:6:737)), [-1:-6: -737])


figure

plot(saltt, [-1:-1: -737]) , hold on
plot(smooth(saltt(1:6:737)), [-1:-6: -737])



