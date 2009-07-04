%%% Reads ouput of Jon's ZEUS code
function [X,Y,Z,RHO,MAGX,MAGY,MAGZ,CURX,CURY,CURZ,CUR] = read_emfdump(fname1,nx,ny,nz);

% Only use the twice inside points to avoid surface defects
nxl = nx-2;
nyl = ny-2;
nzl = nz-2;


% Open file
in = fopen(fname1,'r');

% Read in data
dataemf = fscanf(in,'%e',[19,inf]);
fclose(in);

% Generate Mesgrids
[X,Y,Z] = meshgrid( 1:nx, 1:ny, 1:nz );
RHO = 0.*X;
MAGX = 0.*X;
MAGY = 0.*X;
MAGZ = 0.*X;
CURX = 0.*X;
CURY = 0.*X;
CURZ = 0.*X;
CUR = 0.*X;


for i=1:nxl,
    for j=1:nyl,
        for k=1:nzl,
            ijk = i+1 + nx*( (j+1-1) + ny*(k+1-1) );
            RHO(j,i,k) = dataemf(4,ijk);
            MAGX(j,i,k) = dataemf(14,ijk);
            MAGY(j,i,k) = dataemf(15,ijk);
            MAGZ(j,i,k) = dataemf(16,ijk);
            CURX(j,i,k) = dataemf(17,ijk);
            CURY(j,i,k) = dataemf(18,ijk);
            CURZ(j,i,k) = dataemf(19,ijk);
            CUR(j,i,k) = sqrt(CURX(j,i,k)^2 + CURY(j,i,k)^2 + CURZ(j,i,k)^2);
        end;
    end;
end;
