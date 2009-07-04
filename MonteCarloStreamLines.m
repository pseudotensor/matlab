% Create a Monte-Carlo density distribution of integral curves
%   vector field is (VX,VY,VZ) defined at positions (X,Y,Z) and
%   the streamlines are to be assigned with weighting D
function [hup,hdown] = MonteCarloStreamLines(nlines,X,Y,Z,VX,VY,VZ,D)

%%% Steamline options (up front for convenient changing!)
stp_size = 0.03;
nstps = 500;




sx = zeros(1,nlines);
sy = sx;
sz = sx;
xmax = max(X(:));
xmin = min(X(:));
ymax = max(Y(:));
ymin = min(Y(:));
zmax = max(Z(:));
zmin = min(Z(:));
for i=1:nlines,
    rmc = 1.0;
    Dmc = 0.0;
    while (rmc>Dmc)
        sx(i) = xmin + (xmax-xmin)*rand(1);
        sy(i) = ymin + (ymax-ymin)*rand(1);
        sz(i) = zmin + (zmax-zmin)*rand(1);
        Dmc = interp3(X,Y,Z,D,sx(i),sy(i),sz(i));
        rmc = rand(1);
    end;
end;

hup = streamline(X,Y,Z,VX,VY,VZ,sx,sy,sz,[stp_size, nstps]);
hdown = streamline(X,Y,Z,-VX,-VY,-VZ,sx,sy,sz,[stp_size, nstps]);
