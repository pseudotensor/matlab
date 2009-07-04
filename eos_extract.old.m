function eos_extract()
% nx
% ny
% nc=number of columns
%fid=fopen('/home/jondata/data.head');
%dir='/home/jondata/grmhd-a.9375-rout400new/';
%dir='f:\\matlabscripts\\matlabscripts\\';
%dir='C:\\Documents and Settings\\jon\\My Documents\\';
%dir='C:\\Documents and Settings\\jon\\My Documents\\eoslarge\\';
dir='./';
%dir='C:\\Documents and Settings\\jon\\My Documents\\eossmall\\';
prefix='eos';
prefixo='eosother';
prefix2='eosnew';
suf1='.head';
suf2='.dat';

% SOME PARAMETERS
doclean=1;


% speed of light (cm/s)
c = 2.99792458E10;
% Boltzmann's constant in erg/K
kb = 1.3807E-16;
% Planck's constant
hbar = 1.054592483915517E-27;

numheaders=14;

file1=strcat(dir,prefix,suf1);
file2=strcat(dir,prefix,suf2);
file3=strcat(dir,prefix2,suf2);
file4=strcat(dir,prefixo,suf2);
file5=strcat(dir,prefix2,suf1);

fid=fopen(file1);
[myhead,count]=fscanf(fid,'%g',[numheaders]);
fclose(fid);

ii=1;
nc=myhead(ii); ii=ii+1;
nco=myhead(ii); ii=ii+1;

nrhob=myhead(ii); ii=ii+1;
rhobmin=myhead(ii); ii=ii+1;
rhobmax=myhead(ii); ii=ii+1;

lrhobmin=log10(rhobmin);
lrhobmax=log10(rhobmax);
steplrhob=(lrhobmax-lrhobmin)/(nrhob-1);

ntk=myhead(ii); ii=ii+1;
tkmin=myhead(ii); ii=ii+1;
tkmax=myhead(ii); ii=ii+1;

ltkmin=log10(tkmin);
ltkmax=log10(tkmax);
stepltk=(ltkmax-ltkmin)/(ntk-1);

nhcm=myhead(ii); ii=ii+1;
hcmmin=myhead(ii); ii=ii+1;
hcmmax=myhead(ii); ii=ii+1;

lhcmmin=log10(hcmmin);
lhcmmax=log10(hcmmax);
steplhcm=(lhcmmax-lhcmmin)/(nhcm-1);


ntdyn=myhead(ii); ii=ii+1;
tdynmin=myhead(ii); ii=ii+1;
tdynmax=myhead(ii); ii=ii+1;

ltdynmin=log10(tdynmin);
ltdynmax=log10(tdynmax);
stepltdyn=(ltdynmax-ltdynmin)/(ntdyn-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open eos.dat

fid=fopen(file2);
[mydata,count]=fscanf(fid,'%g',[nc*nrhob*ntk*nhcm*ntdyn]);
fclose(fid);

% make vector into 5-D array
temp=reshape(mydata,[nc,nrhob,ntk,nhcm,ntdyn]);
clear mydata;
% set primary matricies to be each column(field) read in from SM data (i.e.
% read data is setup with column / nx / ny / nz order
mynewdata=permute(temp,[2,3,4,5,1]);
clear temp;
%
% so now all functions are F(rhob,tk,hcm,tdyn)
%

% {tdyn, hcm, rhob, tk, etae, npratio, p_tot, rho_tot, s_tot, Qm, p_photon,
% p_eleposi,p_N,p_nu}
ii=1;
rhob = mynewdata(:,:,:,:,ii); ii=ii+1;
tk = mynewdata(:,:,:,:,ii); ii=ii+1;
hcm = mynewdata(:,:,:,:,ii); ii=ii+1;
tdyn = mynewdata(:,:,:,:,ii); ii=ii+1;
ptot = mynewdata(:,:,:,:,ii); ii=ii+1;
utot = mynewdata(:,:,:,:,ii); ii=ii+1;
stot = mynewdata(:,:,:,:,ii); ii=ii+1;
qm = mynewdata(:,:,:,:,ii); ii=ii+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% some computed things
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total internal energy
%utot = rhotot - rhob;
chi = utot + ptot;
rhobcsq = rhob.*c.*c;
% hspec is dimensionless
hspec = (rhobcsq+chi)./rhobcsq;




% stot is entropy density (erg/K/cc)
% so compute sound speed using specific entropy:
% below makes no sense with pure photons
sspec = stot./(rhobcsq);
% can do below, but ideal gas case gets complicated for U[rho0,Sden]
%sspec = stot./(rhobcsq+utot);

%x=log10(rhob(:,:,:,8));
%y=log10(tk(:,:,:,8));
%z=log10(npratio(:,:,:,8)+1);
%figure; contour(x,y,z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open q_volume.dat

if 0

fid=fopen(file4);
[mydata,count]=fscanf(fid,'%g',[nco*nrhob*ntk*nhcm*ntdyn]);
fclose(fid);

% make vector into 5-D array
temp=reshape(mydata,[nco,ntk,nrhob,nhcm,ntdyn]);
clear mydata;
% set primary matricies to be each column(field) read in from SM data (i.e.
% read data is setup with column / nx / ny / nz order
mynewdata=permute(temp,[2,3,4,5,1]);
clear temp;
%
%
% so now all functions are F(rhob,tk,hcm,tdyn)
%

%


% {tdyn, hcm, rhob, tk, etae, tautel,tautmu,tauael,tauamu,
% qminusel,qminusmu}
ii=1;
etae = mynewdata(:,:,:,:,ii); ii=ii+1;
xnuc = mynewdata(:,:,:,:,ii); ii=ii+1;
npratio = mynewdata(:,:,:,:,ii); ii=ii+1;

pphoton = mynewdata(:,:,:,:,ii); ii=ii+1;
peleposi = mynewdata(:,:,:,:,ii); ii=ii+1;
pN = mynewdata(:,:,:,:,ii); ii=ii+1;
pnu = mynewdata(:,:,:,:,ii); ii=ii+1;

rhophoton = mynewdata(:,:,:,:,ii); ii=ii+1;
rhoeleposi = mynewdata(:,:,:,:,ii); ii=ii+1;
rhoN = mynewdata(:,:,:,:,ii); ii=ii+1;
rhonu = mynewdata(:,:,:,:,ii); ii=ii+1;

sphoton = mynewdata(:,:,:,:,ii); ii=ii+1;
seleposi = mynewdata(:,:,:,:,ii); ii=ii+1;
sN = mynewdata(:,:,:,:,ii); ii=ii+1;
snu = mynewdata(:,:,:,:,ii); ii=ii+1;

tauael = mynewdata(:,:,:,:,ii); ii=ii+1;
taus = mynewdata(:,:,:,:,ii); ii=ii+1;
tautel = mynewdata(:,:,:,:,ii); ii=ii+1;
tautmu = mynewdata(:,:,:,:,ii); ii=ii+1;
tauamu = mynewdata(:,:,:,:,ii); ii=ii+1;

Qmel = mynewdata(:,:,:,:,ii); ii=ii+1;
Qmmu = mynewdata(:,:,:,:,ii); ii=ii+1;
Qmtau = mynewdata(:,:,:,:,ii); ii=ii+1;

qminusel = mynewdata(:,:,:,:,ii); ii=ii+1;
qminusmu = mynewdata(:,:,:,:,ii); ii=ii+1;

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% First construct 1-D grid of new dependent quantities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factor=10;
%factor=1;


% grid of internal energy (utot), pressure (ptot), and \chi=u+p
nu=ntk*factor;
%lumin=log10(rhobmin*c*c);
%lumax=log10(rhobmax*c*c);
lumin=log10(min(min(min(min(utot)))));
lumax=log10(max(max(max(max(utot)))));
lchimin=log10(min(min(min(min(chi)))));
lchimax=log10(max(max(max(max(chi)))));
lpmin=log10(min(min(min(min(ptot)))));
lpmax=log10(max(max(max(max(ptot)))));

lumin=min(min(lumin,lchimin),lpmin);
lumax=max(max(lumax,lchimax),lpmax);

steplu = (lumax-lumin)/(nu-1);
lugrid=lumin:steplu:lumax;
ugrid = 10.^lugrid;

% grid of \chi=u+p
%nchi=ntk*factor;
%lchimin=log10(min(min(min(min(chi)))));
%lchimax=log10(max(max(max(max(chi)))));
%steplchi = (lchimax-lchimin)/(nchi-1);
%lchigrid=lchimin:steplchi:lchimax;
%chigrid = 10.^lchigrid;

% entropy as a variable is only used temporarily
nsspec=ntk*factor;
lsspecmin=log10(min(min(min(min(sspec)))));
lsspecmax=log10(max(max(max(max(sspec)))));
steplsspec = (lsspecmax-lsspecmin)/(nsspec-1);
lsspecgrid=lsspecmin:steplsspec:lsspecmax;
sspecgrid = 10.^lsspecgrid;

% OUTPUT grid of functions of : 1) internal energy (utot) , 2) pressure (ptot), and 3) \chi=u+p
nuout=ntk;
%lumin=log10(rhobmin*c*c);
%lumax=log10(rhobmax*c*c);
luoutmin=lumin;
luoutmax=lumax;
stepluout = (luoutmax-luoutmin)/(nuout-1);
luoutgrid=luoutmin:stepluout:luoutmax;
uoutgrid = 10.^luoutgrid;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now interpolate to obtain desired functions
%
% interpolate in log10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% say X, Y, Xi and get Yi
for p=1:nrhob
    for q=1:nhcm
        for r=1:ntdyn

    %           p
    %        q
    %        r
            
            % check for duplicates in utot
%[xc,yc] = consolidator(utot(p,:,q,r),[],'count');
%myone=ones(length(yc),1);
%mydiffutot=sum(abs(abs(yc)-myone))

%[xc,yc] = consolidator(ptot(p,:,q,r),[],'count');
%myone=ones(length(yc),1);
%mydiffptot=sum(abs(abs(yc)-myone))

%[xc,yc] = consolidator(chi(p,:,q,r),[],'count');
%myone=ones(length(yc),1);
%mydiffchi=sum(abs(abs(yc)-myone))

%[xc,yc] = consolidator(stot(p,:,q,r),[],'count',.05);
%myone=ones(length(yc),1);
%mydiffstot=sum(abs(abs(yc)-myone))

%[xc,yc] = consolidator(log10(stot(p,:,q,r)),[],'count',.05);
%myone=ones(length(yc),1);
%mydiffstot=sum(abs(abs(yc)-myone))


%fprintf('%21.15g\n',log10(stot(p,:,q,r)))



% make sure things interpolating are unique

[lutotptotx,lutotptoty] = consolidator(log10(utot(p,:,q,r)),log10(ptot(p,:,q,r)),'mean',1E-14);

[lhspecx,lhspecy] = consolidator(log10(utot(p,:,q,r)),log10(hspec(p,:,q,r)),'mean',1E-14);

[lptotutotx,lptotutoty] = consolidator(log10(ptot(p,:,q,r)),log10(utot(p,:,q,r)),'mean',1E-14);

[lchiptotx,lchiptoty] = consolidator(log10(chi(p,:,q,r)),log10(ptot(p,:,q,r)),'mean',1E-14);

[lsspecptotx,lsspecptoty] = consolidator(log10(sspec(p,:,q,r)),log10(ptot(p,:,q,r)),'mean',1E-14);
[lsspecutotx,lsspecutoty] = consolidator(log10(sspec(p,:,q,r)),log10(utot(p,:,q,r)),'mean',1E-14);

[lutotqmx,lutotqmy] = consolidator(log10(utot(p,:,q,r)),log10(qm(p,:,q,r)),'mean',1E-14);

[lutotstotx,lutotstoty] = consolidator(log10(utot(p,:,q,r)),log10(stot(p,:,q,r)),'mean',1E-14);
            
%  f(x) = interp(x(t),f(t),x_i)
        
% below is P(rho0,u,H)
% below directly used as P(rho0,u)
% dP/du |rho0
% dP/drho0 | u
PofU(p,:,q,r) = 10.^(myinterp1(1,lutotptotx, lutotptoty, lugrid','linear'));

HofU(p,:,q,r) = 10.^(myinterp1(2,lhspecx, lhspecy, lugrid','linear'));

% below is u(rho0,P,H)
% below directly used as u(rho0,p)
UofP(p,:,q,r) = 10.^(myinterp1(3,lptotutotx, lptotutoty, lugrid','linear'));

% below is P(rho0,\chi,,H)
% below used for P[rho0,chi]
% 1/ (dchi/dp)|rho0
% 1/(drho0/dp)|chi
PofCHI(p,:,q,r) = 10.^(myinterp1(4,lchiptotx, lchiptoty, lugrid','linear'));

% Below is PofS(rho0,S,H)
% below used for c_s^2 = 1/h dp/drho0|S
PofS(p,:,q,r) = 10.^(myinterp1(5,lsspecptotx, lsspecptoty, lsspecgrid','linear'));
UofS(p,:,q,r) = 10.^(myinterp1(6,lsspecptotx, lsspecutoty, lsspecgrid','linear'));

% Below is qmofU(rho0,u,H)
qmofU(p,:,q,r) = 10.^(myinterp1(7,lutotqmx, lutotqmy, lugrid','linear'));

% Below is SofU(rho0,u,H)
SofU(p,:,q,r) = 10.^(myinterp1(8,lutotstotx, lutotstoty, lugrid','linear'));

        end
    end
end

%HofU(1,14,1,1)
%HofU(1,:,1,1)

if doclean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Adjust quantities to be physical when interpolated range is beyond
% existing range
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rhobgrid, ugrid, and sspecgrid need to be turned into full 4-D quantity for cleanvar

% say X, Y, Xi and get Yi

for p=1:nu  % same # when second variable is either U, P, CHI, or S
    for q=1:nhcm
        for r=1:ntdyn
            rhobgrid4d(:,p,q,r) = rhob(:,1,1,1);
            rhobcsqgrid4d(:,p,q,r) = rhobcsq(:,1,1,1);
        end
    end
end
for p=1:nrhob
    for q=1:nhcm
        for r=1:ntdyn
           ugrid4d(p,:,q,r) = ugrid(:);
           sspecgrid4d(p,:,q,r) = sspecgrid(:);
        end
    end
end

PofU = cleanvar(1, PofU, rhobgrid4d, ugrid4d);
HofU = cleanvar(2, HofU, rhobcsqgrid4d, ugrid4d); % specific enthalpy is dimensionless
UofP = cleanvar(3, UofP, rhobgrid4d, ugrid4d);
PofCHI = cleanvar(4, PofCHI, rhobgrid4d, ugrid4d);
PofS = cleanvar(5, PofS, rhobcsqgrid4d, sspecgrid4d); % here S is specific entropy and uses rhobcsq
UofS = cleanvar(6, UofS, rhobcsqgrid4d, sspecgrid4d); % here S is specific entropy and uses rhobcsq
qmofU = cleanvar(7, qmofU, rhobgrid4d, ugrid4d);
SofU = cleanvar(8, SofU, rhobcsqgrid4d, ugrid4d); % here S is entropy density and uses rhobcsq

%totalnan=sum(sum(sum(sum(isnan(UofS)))))
end

%HofU(1,14,1,1)
%HofU(1,:,1,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now compute derivatives using interpolated functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first get rho(i) and S(i)
% assumes rhob is same for all other indicies (true)
 rhocsqi = c.*c.*rhob(:,1,1,1)';
 drhocsqi = gradient(rhocsqi);

 rhoi = rhob(:,1,1,1)';
 drhoi = gradient(rhoi);

 % entropy density (stot) as independent quantity is setup to be same as ugrid
% stoti = ugrid;
% dstoti=gradient(stoti);

 % specific entropy as independent quantity
% sspeci = ugrid./rhocsqi;
 sspeci = sspecgrid;
 dsspeci=gradient(sspeci);

 utoti = ugrid;
 dutoti=gradient(utoti);
 
 chii = ugrid;
 dchii=gradient(chii);
 
% say X, Y, Xi and get Yi
for q=1:nhcm
    for r=1:ntdyn
% NOTE
% The first output FX is always the gradient along the 2nd dimension of F, going across columns. The second output FY is always the gradient along the 1st dimension of F, going across rows. For the third output FZ and the outputs that follow, the Nth output is the gradient along the Nth dimension of F. 

% if PofS(rho0,S,H) then derivative is dS, drho, dH
% gradient gives df/di along each direction
% so if want dP/dS = dP/di / dS/di

% note order of dsspeci and drhocsqi is such that corresponds to output
% derivatives that are flipped for 1 and 2 (row/col) as described above

% when vector is given for second+ terms in gradient, assumed to be
% positions instead of spacing as when scalar input

 [dPofSdS(:,:,q,r), dPofSdrho0(:,:,q,r)] = gradient(PofS(:,:,q,r),sspeci,rhocsqi);

 % PofU(rho0,U,H) then derivative is dU, drho0, dH
 [dPofUdu(:,:,q,r), dPofUdrho0(:,:,q,r)] = gradient(PofU(:,:,q,r),utoti,rhocsqi);

  % PofCHI(rho0, chi ,H) then derivative is dchi, drho0, dH
 [dPofCHIdchi(:,:,q,r), dPofCHIdrho0(:,:,q,r)] = gradient(PofCHI(:,:,q,r),chii,rhocsqi);

 % SofU(rho0,U,H) then derivative is dU, drho0, dH
 [dSofUdu(:,:,q,r), dSofUdrho0(:,:,q,r)] = gradient(SofU(:,:,q,r),utoti,rhocsqi);
    end
end

%gradPofCHI=gradient(ptot(1,:,1,1),chi(1,:,1,1));
%figure; plot(chi(1,:,1,1),gradPofCHI)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Some derivatives need to be have a change of variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% say X, Y, Xi and get Yi
for p=1:nrhob
    for q=1:nhcm
        for r=1:ntdyn
            
% below is presently dP/drho(rho0,S)|S and need dP/drho(rho0,u)|S so change S->u
% dPofSdrho0(:,:,q,r)

% GODMARK: Should below be utot(rho0,T) -> UofS(rho0,S) ?  Seems so.

[lutotdpx,lutotdpy] = consolidator(log10(UofS(p,:,q,r)),dPofSdrho0(p,:,q,r),'mean',1E-14);

% Below is dP/drho0(U,rho0,H)|S
% below used for c_s^2 = 1/h dp/drho0|S
dPofUdrho0cS(p,:,q,r) = myinterp1(9,lutotdpx, lutotdpy, lugrid','linear');

        end
    end
end


if doclean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Need to clean dPofUdrho0cS since reinterpolated from F(rho0,S) back to F(rho0,U)
%
% Notice that this uses rhobcsq instead of rhob because derivatives are
% assumed to in dimensionless ratios
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dPofUdrho0cS = cleanvar(9, dPofUdrho0cS, rhobcsqgrid4d, ugrid4d);

end


            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now compute things that don't involve derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c_s^2 in (cm/s)^2
% dPofUdrho0cS already has 1/c^2 in it
% so below is per unit c^2
cs2 = (1.0./HofU).*(dPofUdrho0cS);
cs2cgs = cs2.*c.*c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% list of things to output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PofU(rho0,U)
% UofP(rho0,P)

% dPofUdrho0(rho0,U,H)
% dPofUdu(rho0,U,H)

% cs2(rho0,U,H)

% sofU(rho0,U,H)
% dSofUdrho0(rho0,U,H)
% dSofUdu(rho0,U,H)

% PofCHI(rho0,CHI,H)
% dPofCHIdrho0(rho0,CHI,H)
% dPofCHIdchi(rho0,CHI,H)


% qmofU(rho0,u,H)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Downsample for writing file
%
% downsample: ugrid since all F(:,U,:,:) are of larger-than-desired size
%
% for quantities that have large dynamic range use log interpolation
% unless quantities can naturally be 0 or negative such as derivatives or
% things from derivatives (e.g. cs2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:nrhob
    for q=1:nhcm
        for r=1:ntdyn
   
%            [lutotdpx,lutotdpy] = consolidator(log10(UofS(p,:,q,r)),log10(dPofSdrho0(p,:,q,r)),'mean',1E-14);

 PofUout(p,:,q,r) = 10.^(myinterp1(1,lugrid, log10(PofU(p,:,q,r)), luoutgrid','linear'));
 HofUout(p,:,q,r) = 10.^(myinterp1(2,lugrid, log10(HofU(p,:,q,r)), luoutgrid','linear'));
 UofPout(p,:,q,r) = 10.^(myinterp1(3,lugrid, log10(UofP(p,:,q,r)), luoutgrid','linear'));
% dPofUdrho0out(p,:,q,r) = 10.^(myinterp1(20,lugrid, log10(dPofUdrho0(p,:,q,r)), luoutgrid','linear'));
% dPofUduout(p,:,q,r) = 10.^(myinterp1(21,lugrid, log10(dPofUdu(p,:,q,r)), luoutgrid','linear'));
 dPofUdrho0out(p,:,q,r) = myinterp1(20,lugrid, dPofUdrho0(p,:,q,r), luoutgrid','linear');
 dPofUduout(p,:,q,r) = myinterp1(21,lugrid, dPofUdu(p,:,q,r), luoutgrid','linear');
 cs2out(p,:,q,r) = myinterp1(22,lugrid, cs2(p,:,q,r), luoutgrid','cubic');
 cs2cgsout(p,:,q,r) = myinterp1(23,lugrid, cs2cgs(p,:,q,r), luoutgrid','cubic');
 SofUout(p,:,q,r) = 10.^(myinterp1(8,lugrid, log10(SofU(p,:,q,r)), luoutgrid','linear'));
% dSofUdrho0out(p,:,q,r) = 10.^(myinterp1(24,lugrid, log10(dSofUdrho0(p,:,q,r)), luoutgrid','linear'));
% dSofUduout(p,:,q,r) = 10.^(myinterp1(25,lugrid, log10(dSofUdu(p,:,q,r)), luoutgrid','linear'));
 dSofUdrho0out(p,:,q,r) = myinterp1(24,lugrid, dSofUdrho0(p,:,q,r), luoutgrid','linear');
 dSofUduout(p,:,q,r) = myinterp1(25,lugrid, dSofUdu(p,:,q,r), luoutgrid','linear');
 PofCHIout(p,:,q,r) = 10.^(myinterp1(4,lugrid, log10(PofCHI(p,:,q,r)), luoutgrid','linear'));
% dPofCHIdrho0out(p,:,q,r) = 10.^(myinterp1(26,lugrid, log10(dPofCHIdrho0(p,:,q,r)), luoutgrid','linear'));
% dPofCHIdchiout(p,:,q,r) = 10.^(myinterp1(27,lugrid, log10(dPofCHIdchi(p,:,q,r)), luoutgrid','linear'));
 dPofCHIdrho0out(p,:,q,r) = myinterp1(26,lugrid, dPofCHIdrho0(p,:,q,r), luoutgrid','linear');
 dPofCHIdchiout(p,:,q,r) = myinterp1(27,lugrid, dPofCHIdchi(p,:,q,r), luoutgrid','linear');
  
 qmofUout(p,:,q,r) = 10.^(myinterp1(28,lugrid, log10(qmofU(p,:,q,r)), luoutgrid','linear'));
 
        end
    end
end

%HofUout(1,:,1,1)



clear HofU cs2 PofU UofP dPofUdrho0 dPofUdu cs2cgs SofU dSofUdrho0 dSofUdu PofCHI dPofCHIdrho0 dPofCHIdchi qmofU;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Adjust quantities to be physical in case of numerical error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[mincs2out,mincs2outI]=min(min(min(min(cs2out))))
[maxcs2out,maxcs2outI]=max(max(max(max(cs2out))))

[mincs2cgsout,mincs2cgsoutI]=min(min(min(min(cs2cgsout))))
[maxcs2cgsout,maxcs2cgsoutI]=max(max(max(max(cs2cgsout))))

% adjust speed of sound to be no smaller than 0 and no larger than 1

for p=1:ntdyn
    for o=1:nhcm
        for n=1:nuout
            for m=1:nrhob
                    
                    if cs2out(m,n,o,p)>1.0
                       cs2out(m,n,o,p)=1.0-1E-14;
                    end
                    if cs2out(m,n,o,p)<0.0
                       cs2out(m,n,o,p)=0.0;
                    end

                    if cs2cgsout(m,n,o,p)>c.*c
                       cs2cgsout(m,n,o,p)=(1.0-1E-14).*c.*c;
                    end
                    if cs2cgsout(m,n,o,p)<0.0
                       cs2cgsout(m,n,o,p)=0.0;
                    end

            end
        end
    end
end

[mincs2out,mincs2outI]=min(min(min(min(cs2out))))
[maxcs2out,maxcs2outI]=max(max(max(max(cs2out))))

[mincs2cgsout,mincs2cgsoutI]=min(min(min(min(cs2cgsout))))
[maxcs2cgsout,maxcs2cgsoutI]=max(max(max(max(cs2cgsout))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% write data to file
%
% Note that all ratios (or derivatives) are dimensionless
%
% all other quantities, including sound speed, have dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% GODMARK: Check units of entropy stuff

%                HofUout(m,n,o,p) , cs2out(m,n,o,p), ...
% %21.15g %21.15g 


%%%%%%%%%%%%%%%%%%%%%%%%
%
% Notice that below indicies are in C-order and C-style (start with 0
% instead of 1)
%
%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen(file3,'w');
for p=1:ntdyn
    for o=1:nhcm
        for n=1:nuout
            for m=1:nrhob
                fprintf(fid,'%3d %3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g ', ...
                m-1, n-1, o-1, p-1, ...
                rhob(m,n,o,p), uoutgrid(n), hcm(m,n,o,p), tdyn(m,n,o,p), ...
                PofUout(m,n,o,p), UofPout(m,n,o,p), ...
                dPofUdrho0out(m,n,o,p), dPofUduout(m,n,o,p), ...
                cs2cgsout(m,n,o,p), ...
                SofUout(m,n,o,p), dSofUdrho0out(m,n,o,p), dSofUduout(m,n,o,p), ...
                PofCHIout(m,n,o,p), dPofCHIdrho0out(m,n,o,p), dPofCHIdchiout(m,n,o,p), ...
                qmofUout(m,n,o,p));
                fprintf(fid,'\n');
            end
        end
    end
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% write header to file
%
% Note that steps and all such things are computed same as in Kaz's code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fid=fopen(file5,'w');
fprintf(fid,'%d\n',20); % number of output colums
fprintf(fid,'%d %21.15g %21.15g %21.15g\n',nrhob,lrhobmin,lrhobmax,steplrhob);
%fprintf(fid,'%d %21.15g %21.15g %21.15g\n',ntk,ltkmin,ltkmax,stepltk);
fprintf(fid,'%d %21.15g %21.15g %21.15g\n',nuout,luoutmin,luoutmax,stepluout);
fprintf(fid,'%d %21.15g %21.15g %21.15g\n',nhcm,lhcmmin,lhcmmax,steplhcm);
fprintf(fid,'%d %21.15g %21.15g %21.15g\n',ntdyn,ltdynmin,ltdynmax,stepltdyn);
fclose(fid);


fclose('all');


quit;


end


% 1:  log10 PofU (rho0,U)
% 2:  log10 HofU (rho0,U)
% 3:  log10 UofP (rho0,P)
% 4:  log10 PofCHI (rho0,CHI)
% 5:  log10 PofS (rho0,S)               % S is specific entropy
% 6:  log10 UofS (rho0,S)               % S is specific entropy
% 7:  log10 qmofU (rho0,U)
% 8:  log10 SofU (rho0,U)               % S is entropy density
% 9:  lin   dPofUdrho0cS (rho0,U) | S   % S is specific entropy
% 20: lin   dPofUdrho0out (rho0,U)
% 21: lin   dPofUduout (rho0,U)
% 22: lin   cs2/csq (rho0,U)
% 23: lin   cs2cgs (rho0,U)
% 24: lin   dSofUdrho0out (rho0,U)  % S is entropy density
% 25: lin   dSofUduout (rho0,U)     % S is entropy density
% 26: lin   dPofCHIdrho0out (rho0,CHI)
% 27: lin   dPofCHIdchiout (rho0,CHI)
% 28: log10 qmofUout (rho0,U)

% input is full 4-D quantities
function funout = cleanvar(whichvar, funin, X, Y)

myisnan=isnan(funin);

if whichvar==1
  funin(myisnan) = idealPofU(X(myisnan),Y(myisnan));
end
if whichvar==2
  funin(myisnan) = idealHofU(X(myisnan),Y(myisnan));
end
if whichvar==3
  funin(myisnan) = idealUofP(X(myisnan),Y(myisnan));
end
if whichvar==4
  funin(myisnan) = idealPofCHI(X(myisnan),Y(myisnan));
end
if whichvar==5
  funin(myisnan) = idealPofS(X(myisnan),Y(myisnan));
end
if whichvar==6
  funin(myisnan) = idealUofS(X(myisnan),Y(myisnan));
end
if whichvar==7
  funin(myisnan) = idealqmofU(X(myisnan),Y(myisnan));
end
if whichvar==8
  funin(myisnan) = idealSdenofU(X(myisnan),Y(myisnan));
end
if whichvar==9
   funin(myisnan) = idealdPofUdrho0cS(X(myisnan),Y(myisnan));
end

% assume rest are derivative quantities (i.e. derived from the above)


% assume we fixed funin and set to funout
funout = funin;

end



function varargout = myinterp1(whichvar, varargin)


%extrap1(lsspecptotx, lsspecptoty, lsspecgrid','linear');

% do extrapolation
% as designed, derived quantities may obtain awkward values

%if 0


%if whichvar==8

% NORMAL:
if whichvar==5 || whichvar==6

% DEBUG:
%if whichvar==5 || whichvar==6 || whichvar==1 || whichvar==2 || whichvar==3 || whichvar==4 || whichvar==7 || whichvar==8 || ...
%        whichvar == 20 || whichvar == 21 || whichvar==22 || whichvar==23 || whichvar==24 || whichvar==25 || whichvar==26 || whichvar==27 || whichvar==28

% DEBUG:
% if 1 || whichvar==5 || whichvar==6
    varargout{:}=extrap1(varargin{:});
else
    
    %whichvar
    
% don't do extrapolation
% returns NaN's so can extrapolate PER quantity rather than deriving things
% from extrapolated quantities
    varargout{:}=interp1(varargin{:});
end
% try just setting NaN's to 0
%varargout{1}(isnan(varargout{1}(:))) = 0;

end



function out = idealPofU(rho0, U)

GAMMA=(4.0/3.0);

out = (GAMMA-1.0).*U;

end

% ratios of energies assumed to be dimensionless
function out = idealHofU(rhobcsq, U)

GAMMA=(4.0/3.0);
out = 1.0 + GAMMA.*U./(rhobcsq);

end


function out = idealUofP(rho0, P)

GAMMA=(4.0/3.0);

out = P./(GAMMA-1.0);

end

function out = idealPofCHI(rho0, CHI)

GAMMA=(4.0/3.0);

out = ((GAMMA-1.0)./GAMMA).*CHI;

end



% Sden 
function out = idealPofSden(rhobcsq, Sden)

GAMMA=(4.0/3.0);
C = 0.00039324313707644323;
%out = (rho0.^GAMMA) .* exp(Sden./(kb.*rho0./mN));

out = (GAMMA-1.0).*((Sden./C).^(GAMMA));


end



% Sden
function out = idealUofSden(rhobcsq, Sden)

GAMMA=(4.0/3.0);

out = idealPofSden(rhobcsq,Sden)./(GAMMA-1.0);

%out = (Sden./C).^(4.0/3.0);

end


% here S is specific entropy
function out = idealPofS(rhobcsq, S)

% consistent with how defined specific entropy above (sspec)
Sden = S.*rhobcsq;

out = idealPofSden(rhobcsq, Sden);

end



% here S is specific entropy
function out = idealUofS(rhobcsq, S)

GAMMA=(4.0/3.0);

% idealPofS takes in specific entropy
out = idealPofS(rhobcsq,S)./(GAMMA-1.0);

end





function out = idealqmofU(rho0, U)

GAMMA=(4.0/3.0);

% put in something small so can use log-log interpolation
out = 1E-20;

end



% S = entropy density
function out = idealSdenofU(rhobcsq, U)

GAMMA=(4.0/3.0);
C = 0.00039324313707644323;

% from Kaz's formula for photon entropy and see kaz_entropy_conversion.nb
% for how converted to cgs
out = C.*U.^(3.0/4.0);

end



% S = specific entropy
function out = idealSofU(rhobcsq, U)

GAMMA=(4.0/3.0);
P = (GAMMA-1.0)*U;
Sden = idealSdenofU(rhobcsq, U);

% true specific entropy?
%out = Sden./(rhobcsq + U + P);
% makes no sense if no baryons
out = Sden./rhobcsq;

end



% ratios from derivatives assumed to be dimensionless
% holding constant: S = specific entropy
function out = idealdPofUdrho0cS(rhobcsq,U)

GAMMA=(4.0/3.0);
hspec = idealHofU(rhobcsq,U);
P = (GAMMA-1.0).*U;

% dP/drho|S = h c_s^2 = h GAMMA P/(\rho_0+u+P) = h GAMMA(GAMMA-1)u/(\rho_0+u+P)
out = hspec.*GAMMA.*P./(rhobcsq+U+P);

end




