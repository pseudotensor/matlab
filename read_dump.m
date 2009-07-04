%%% Reads ouput of Jon's ZEUS code
function [X,Y,Z,RHO,U,PHI,VX,VY,VZ,BX,BY,BZ] = read_dump(fname)

% Open file
in=fopen(fname,'r');

% Read in header information
fgetl(in);
data = fscanf(in,'%e',[3,1]);
fgetl(in);
data = fscanf(in,'%e',[4,1]);
fgetl(in);
data = fscanf(in,'%i',[4,1]);
nx = data(1);
ny = data(2);
nz = data(3);
fgetl(in);

% Read in data
data = fscanf(in,'%e',[9,inf]);
fclose(in);

% Generate Mesgrids
[X,Y,Z] = meshgrid( 1:nx, 1:ny, 1:nz );
RHO = 0.*X;
U = 0.*X;
PHI = 0.*X;
VX = 0.*X;
VY = 0.*X;
VZ = 0.*X;
BX = 0.*X;
BY = 0.*X;
BZ = 0.*X;

for i=1:nx,
    for j=1:ny,
        for k=1:nz,
            ijk = i + nx*( (j-1) + ny*(k-1) );
            RHO(i,j,k) = data(1,ijk);
            U(i,j,k) = data(2,ijk);
            PHI(i,j,k) = data(3,ijk);
            VX(i,j,k) = data(4,ijk);
            VY(i,j,k) = data(5,ijk);
            VZ(i,j,k) = data(6,ijk);
            BX(j,i,k) = data(7,ijk);
            BY(j,i,k) = data(8,ijk);
            BZ(j,i,k) = data(9,ijk);
        end;
    end;
end;
