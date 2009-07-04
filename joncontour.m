function [C1]=joncontour(nx,ny)
% nx=number of columns
% ny=number of rows
fid=fopen('c:\Documents and Settings\jon.JONXP\My Documents\research\r.txt');
[r,count]=fscanf(fid,'%f',[nx,ny]);
fclose(fid);
fid=fopen('c:\Documents and Settings\jon.JONXP\My Documents\research\h.txt');
[h,count]=fscanf(fid,'%f',[nx,ny]);
fclose(fid);
fid=fopen('c:\Documents and Settings\jon.JONXP\My Documents\research\i.txt');
[i,count]=fscanf(fid,'%f',[nx,ny]);
fclose(fid);
fid=fopen('c:\Documents and Settings\jon.JONXP\My Documents\research\j.txt');
[j,count]=fscanf(fid,'%f',[nx,ny]);
fclose(fid);
fid=fopen('c:\Documents and Settings\jon.JONXP\My Documents\research\aphi.txt');
[aphi,count]=fscanf(fid,'%f',[nx,ny]);
fclose(fid);
fid=fopen('c:\Documents and Settings\jon.JONXP\My Documents\research\w.txt');
[w,count]=fscanf(fid,'%f',[nx,ny]);
fclose(fid);

%r=transpose(r)
%h=transpose(h)
%i=transpose(i)
%j=transpose(j)
%aphi=transpose(aphi)
%w=transpose(w)

%size(aphi)
%figure; contour(r,h,aphi,50); title(['aphi']);
%figure; contour(r,w,aphi,50); title(['aphi']);
figure; [C,h]=contour(r,w,aphi,[0.002,0.002]); title(['aphi']);


%figure; contour(C); title(['C']);

C
[sizex sizey]=size(C(1,:))

sizex
sizey

myr=C(1,2:sizey)
myw=C(2,2:sizey)

mymat=zeros(2, sizey-1);
mymat(1,:)=myr
mymat(2,:)=myw
%C(2)



figure; plot(myr,myw); title(['aphi3']);
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

fid=fopen('C:\Documents and Settings\jon.JONXP\My Documents\research\C1.txt','w');
for m = 1:sizey-1
    for n = 1:2
        fprintf(fid,'%f ',mymat(n,m))
    end
    fprintf(fid,'\n')
end
fclose(fid)

end
