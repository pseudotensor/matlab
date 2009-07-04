function [C1,C2]=joncorr(nx,ny)
% nx=number of columns
% ny=number of rows
fid=fopen('vsaclean.txt');
[A,count]=fscanf(fid,'%f',[nx,ny])
fclose(fid)

B=transpose(A)

C1 = corrcoef(A)
C2 = corrcoef(B)

fid=fopen('C:\Documents and Settings\jon.JONXP\My Documents\research\bztables\C1.txt','w');
for m = 1:ny
    for n = 1:ny
        fprintf(fid,'%f ',C1(n,m))
    end
    fprintf(fid,'\n')
end
fclose(fid)

fid=fopen('C:\Documents and Settings\jon.JONXP\My Documents\research\bztables\C2.txt','w');
for m = 1:nx
    for n = 1:nx
        fprintf(fid,'%f ',C2(n,m))
    end
    fprintf(fid,'\n')
end
fclose(fid)


end
