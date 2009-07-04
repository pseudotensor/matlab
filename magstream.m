%%% Function to generate magnetic field lines from Jon's ZEUS code output
function magstream(nlines,fname)

%% (A) Clear current figure if need be
clf;

%% (B) Generate vector field
[X,Y,Z,RHO,U,PHI,VX,VY,VZ,BX,BY,BZ] = read_dump(fname);

%% (C)
p = patch(isosurface(X,Y,Z,log10(RHO),-3));
isonormals(X,Y,Z,log10(RHO),p);
set(p,'FaceColor','green','EdgeColor','none','FaceAlpha',0.5);
camlight;
lighting gouraud;

Bnorm = sqrt(BX.^2+BY.^2+BZ.^2);
Bnorm = Bnorm/max(Bnorm(:));
[bhup,bhdown] = MonteCarloStreamLines(nlines,X,Y,Z,BX,BY,BZ,Bnorm);
set(bhup,'LineWidth',1,'Color','red');
set(bhdown,'LineWidth',1,'Color','blue');


%%   Right now just sample evenly in middle of the plane
[nx ny nz] = size(X);


axis([1 nx 1 ny 1 nz]);
axis square;

set(gcf,'color','black');
set(gca,'color','black','xcolor','white','ycolor','white','zcolor','white');


% TUBE TRIAL
% htubes = streamtube(X,Y,Z,BX,BY,BZ,sx,sy,sz);
% set(htubes,'facecolor','blue','edgecolor','none');
% camlight;
% lighting gouraud;
