function jonstreamline1(wc)
% nx
% ny
% nz
dir='/home/jondata/nsnew1/run/dumps/'
prefix1='myhead'
prefix2='myiaphi'
prefix3='streamlines'
suf1='.txt'

file1=strcat(dir,prefix1,suf1);
file2=strcat(dir,prefix2,suf1);
file3=strcat(dir,prefix3,suf1);


fid=fopen(file1);
[myhead,count]=fscanf(fid,'%g',[16]);
fclose(fid);

time=myhead(1)
nx=myhead(2)
ny=myhead(3)
nz=myhead(3) % assumed
startx=myhead(4)
starty=myhead(5)
startz=myhead(5) % assumed
dx=myhead(6)
dy=myhead(7)
dz=myhead(7) % assumed

fid=fopen(file2);
[mydata,count]=fscanf(fid,'%g',[nx*ny*nz]);
fclose(fid);

% make vector into 3-D array
mynewdata=reshape(mydata,[nx,ny,nz]);
%mynewdata=permute(temp,[2,3,1]);

% mynewdata(R,Z,Y)



% showing Z-Y or Z-R plane
% bad plane for spherical r.
%lower values are bad somehow, only for 64/65
% caused by planar_interpolation not acting properly on one edge
god=reshape(mynewdata(64,:,:),[ny,nz]);
contour(god,[wc,wc]);

% showing Z-Y or Z-R plane
god=reshape(mynewdata(:,:,64),[nx,ny]);
contour(god,[wc,wc]);

% showing R-Y plane at Z=middle
% problem with y or R axis for aphi
god=reshape(mynewdata(:,64,:),[nx,nz]);
contour(god,[wc,wc]);

figure;
p=patch(isosurface(mynewdata,wc));
isonormals(mynewdata,p);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1]);
isocolors(mynewdata,p)
shading interp
axis(volumebounds(mynewdata))
%axis([30 80 30 80 30 80]);
view(3);
camlight;
lighting phong;

% or
%figure;
%h= contourslice(mynewdata,[1:9],[],[0],linspace(0,0.3,10));
%axis(volumebounds(mynewdata));
%camva(24); camproj perspective;
%campos('auto');
%set(gcf,'Color',[0.5,0.5,0.5],'Renderer','zbuffer')
%set(gca,'Color','black','XColor','white')
%box on

figure;
god=reshape(mynewdata(64,:,:),[nx ny]);
[FX,FY]=gradient(god,200/128);
FXnew=-FY;
FYnew=FX;
[sx,sy]=meshgrid(70,70);
h=streamline(FXnew,FYnew,sx,sy);
%h=streamtube(FXnew,FYnew,sx,sy);
set(h,'Color','red');
view(3);
% now do something with data

% mynewdata(R,Z,Y)

% below gives R vs Y
%god=reshape(mynewdata(:,64,:),[nx ny]);
% below gives Y vs. Z
god=reshape(mynewdata(64,:,:),[nx ny]);
% below gives R vs. Z
%god=reshape(mynewdata(:,:,64),[nx ny]);
[FX,FY]=gradient(god);
streamslice(FX,FY,'cubic');

%fid=fopen(file3,'w');
%for m = 1:sizey-1
%    fprintf(fid,'%21.15g ',myx(m));
%    for n=1:nc
%        fprintf(fid,'%21.15g ',contourdata(m,n));
%    end
%    fprintf(fid,'\n');
%end
%fclose(fid);

end
