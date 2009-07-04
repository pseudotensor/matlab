function [C1]=jonfft
% nx=number of columns
% ny=number of rows
fid=fopen('c:\Documents and Settings\jon.JONXP\My Documents\research\myt.txt');
[t,count]=fscanf(fid,'%f');
fclose(fid);
fid=fopen('c:\Documents and Settings\jon.JONXP\My Documents\research\mydm.txt');
[dm,count]=fscanf(fid,'%f');
fclose(fid);

nx=count

figure;
plot(t,dm);
title('dm vs time')
xlabel('time GM/c^3')

P=fft(dm);

halfnx=round(nx/2);
halfnxp1=halfnx+1;

f = (1/2)*(0:halfnx)/nx;
figure;
plot(f,P(1:halfnxp1))
title('Frequency content of y')
xlabel('frequency (Hz)')


%size(aphi)
%contour(aphi)

%fid=fopen('C:\Documents and Settings\jon.JONXP\My Documents\research\C1.txt','w');
%for m = 1:sizey-1
%    for n = 1:2
%        fprintf(fid,'%f ',mymat(n,m))
%    end
%    fprintf(fid,'\n')
%end
%fclose(fid)

end
