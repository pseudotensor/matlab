function [C]=joncorr(nx,ny)
% nx=number of columns
% ny=number of rows
fid=fopen('vsaclean.txt');
[A,count]=fscanf(fid,'%f',[nx,ny])
fclose(fid)

C = corrcoef(A)

end
