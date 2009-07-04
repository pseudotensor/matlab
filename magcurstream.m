%%% Function to generate magnetic field lines from Jon's ZEUS code output
function magcurstream(nblines,nclines,bfname,cfname)

%% (A) Clear current figure if need be
clf;

%% (B) Generate vector field
[X,Y,Z,RHO,U,PHI,VX,VY,VZ,BX,BY,BZ] = read_dump(bfname);
[nx ny nz] = size(X);

%% (C)
p = patch(isosurface(X,Y,Z,log10(RHO),-2));
isonormals(X,Y,Z,log10(RHO),p);
set(p,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.4);
camlight;
lighting gouraud;

%% (D) Make field lines
Bnorm = sqrt(BX.^2+BY.^2+BZ.^2);
Bnorm = Bnorm/max(Bnorm(:));
Bnorm = Bnorm.^2;
[bhup,bhdown] = MonteCarloStreamLines(nblines,X,Y,Z,BX,BY,BZ,Bnorm);
set(bhup,'LineWidth',1,'Color','red');
set(bhdown,'LineWidth',1,'Color','red');

%% (E) Make current lines
[Xc,Yc,Zc,RHOc,MAGXc,MAGYc,MAGZc,CURXc,CURYc,CURZc,CURc] = read_emfdump(cfname,nx,ny,nz);
Cnorm = CURc/max(CURc(:));
Cnorm = Cnorm.^2;
[chup,chdown] = MonteCarloStreamLines(nclines,Xc,Yc,Zc,CURXc,CURYc,CURZc,Cnorm);
set(chup,'LineWidth',1.5,'Color','green');
set(chdown,'LineWidth',1.5,'Color','green');


axis([1 nx 1 ny 1 nz]);
axis square;

set(gcf,'color','black');
set(gca,'color','black','xcolor','white','ycolor','white','zcolor','white');


% TUBE TRIAL
% htubes = streamtube(X,Y,Z,BX,BY,BZ,sx,sy,sz);
% set(htubes,'facecolor','blue','edgecolor','none');
% camlight;
% lighting gouraud;
