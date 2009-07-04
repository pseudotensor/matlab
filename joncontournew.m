function joncontournew()
% nx
% ny
% nc=number of columns
%fid=fopen('/home/jondata/data.head');
%dir='/home/jondata/grmhd-a.9375-rout400new/'
dir='./'
prefix='caphidata'
prefix2='contoursks'
suf1='.head'
suf2='.txt'

file1=strcat(dir,prefix,suf1);
file2=strcat(dir,prefix,suf2);
file3=strcat(dir,prefix2,suf2);

fid=fopen(file1);
[myhead,count]=fscanf(fid,'%g',[5]);
fclose(fid);

whichaphi=myhead(1)
nc=myhead(2)
nx=myhead(3)
ny=myhead(4)
mycount=myhead(5)

fid=fopen(file2);
[mydata,count]=fscanf(fid,'%g',[nx*ny*nc]);
fclose(fid);

% make vector into 3-D array
temp=reshape(mydata,[nc,nx,ny]);
% set primary matricies to be each column(field) read in from SM data (i.e.
% read data is setup with column / nx / ny order
mynewdata=permute(temp,[2,3,1]);

% {myr myh myti mytj myaphi myw1 myw2 myw3 myuu1 myfastv1p myfastv1m mybsq myB1 myB3}
ii=1;
r = mynewdata(:,:,ii); ii=ii+1;
h = mynewdata(:,:,ii); ii=ii+1;
i = mynewdata(:,:,ii); ii=ii+1;
j = mynewdata(:,:,ii); ii=ii+1;
x1 = mynewdata(:,:,ii); ii=ii+1;
x2 = mynewdata(:,:,ii); ii=ii+1;
aphi = mynewdata(:,:,ii); ii=ii+1; % 8 at end
II = find(aphi>1E10);
% ignore out of bounds values setup in SM
aphi(II)=NaN;
% ii now sets starting position for other things

%size(aphi)
%figure; contour(r,h,aphi,50); title(['aphi']);
%figure; contour(r,w,aphi,50); title(['aphi']);
%figure; [C,ha]=contour(r,w1,aphi,[0.002,0.002]); title(['aphi']);

whichx=r;
[Ctest,ha]=contour(whichx,mynewdata(:,:,ii),aphi,[whichaphi,whichaphi]);
[sizex sizey]=size(Ctest(1,:));
sizex
sizey
myx=fliplr(Ctest(1,2:sizey));
contourdata=zeros(sizey-1,nc);


for m=1:nc
    m
    %figure;
    [Ctest,ha]=contour(whichx,mynewdata(:,:,m),aphi,[whichaphi,whichaphi]);
    %title(['m']);    
    contourdata(:,m)=fliplr(Ctest(2,2:sizey));
end


%mymat=zeros(2, sizey-1);
%mymat(1,:)=myr
%mymat(2,:)=myw
%C(2)



%figure; plot(myr,myw); title(['aphi3']);
%figure; plot(C(2,:)); title(['aphi4']);

%figure; imagesc(C(:,:),[nx ny]); title(['aphi2']);



%figure; contour(r,h,w,50); title(['w']);

%aphi(:)=A(5,:,:)

%h=A(2)
%i=A(3)
%j=A(4)
%w=A(6)


%size(aphi)
%contour(aphi)

fid=fopen(file3,'w');
for m = 1:sizey-1
    fprintf(fid,'%21.15g ',myx(m));
    for n=1:nc
        fprintf(fid,'%21.15g ',contourdata(m,n));
    end
    fprintf(fid,'\n');
end
fclose(fid);

quit;

end
