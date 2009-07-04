function eos_extract()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SOME PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set whether on laptop (interactive) or relativity (matlab script)
  onrelativity=1

  % whether to fix utot if using old Kaz code result where I didn't subtract
  % off rhob c^2
  utotfix = 0

  % whether to use analytical fit or numerical values to set degenerate (offset) values
  utotdegenanalytic=0

  % whether to set "0" for log output as degeneracy fitting formula (to be added by in when inside HARM)
  % trying to get better temperature resolution for the low-temperature domain
  utotdegencut = 1


  % whether to force functions to be monotonic as functions of tk
  forcemonotk = 1

  % whether to "clean" solution if out of bounds
  doclean=1;

  % whether to specify a limit log(u) range (for zoom tables or if user knows that automatically determined range is too large for some reason)
  %specifylurange=1

  % whether to compute sound speed before or after interpolation
  % 0=compute after using entropy (HELM's entropy is messed up, but LS and Shen are fine)
  % 1=compute before using only non-entropy data (BEST)
  % 2=compute both and use smaller of two versions if non-negative and use larger if one is negative (avoids difference errors) (shouldn't use for HELM since HELM's entropy is messed up)
  preinterpsoundspeed=0


  %%%%%%%%%%%% TODO:
  % 1) Need to specify base and log-offset for each variable
  % 2) Then need to output a table for each (rho(X)Hcm(X)Tdynorye) the umin, umax, pmin, pmax, chimin, chimax.  This way I can get high accuracy for T~0Kelvin by starting the table near T~0 and resolving the space near T~0.
  %  ACTUALLY: data is already really there.  I just need to use MATLAB to find min/max of each and make sure to have that as first and last data point for each rho,H,T.  Then in HARM I use 0 and N-1 values as min/max.  Only issue is amount of computation using logs and pow's to process that raw data.
  % Should probably still store it separately, but can be done in HARM instead of MATLAB
  % In matlab each N points should range from the min/max found for EACH rho,H,T instead of globally

  % assume for now offset is 0 and base is 10
  % for now use functional offset


  CONTOL=1E-13;

  OUTBOUNDSVALUE=1E-20;



  % nx
  % ny
  % nc=number of columns
  %fid=fopen('/home/jondata/data.head');
  %dir='/home/jondata/grmhd-a.9375-rout400new/';
  %dir='f:\\matlabscripts\\matlabscripts\\';
  %dir='C:\\Documents and Settings\\jon\\My Documents\\';
  if onrelativity==0
    % BELOW USED ON laptop
    %dir='C:\\Documents and Settings\\jon\\My Documents\\grbjet\\eoslarge\\';
    dir='C:\\Documents and Settings\\jon\\My Documents\\grbjet\\eoslargesimple\\';
  end
  if onrelativity==1
    % BELOW USED ON relativity
    dir='./';
  end
  %dir='C:\\Documents and Settings\\jon\\My Documents\\eossmall\\';
  prefix='eos';
  prefixo='eosother';
  prefix2='eosnew';
  prefix3='eosdegennew';
  prefix4='eosparms'
  prefix5='eosmonodegen'
  suf1='.head';
  suf2='.dat';


  % Thompson et al. (2003) way to merge Helmholtz EOS and LS EOS
  % rho<4E7g/cc use Helmholtz
  % rho>6E7g/cc use LS between quadratically interpolate
  % pressure difference order 1%, entropy order 5-10%.





  % speed of light (cm/s)
  c = 2.99792458E10;
  % Boltzmann's constant in erg/K
  kb = 1.3807E-16;
  % Planck's constant
  % hbar = 1.054592483915517E-27;

  
  file1=strcat(dir,prefix,suf1);
  file2=strcat(dir,prefix,suf2);
  file3=strcat(dir,prefix2,suf2);
  file4=strcat(dir,prefixo,suf2);
  file5=strcat(dir,prefix2,suf1);
  file6=strcat(dir,prefix3,suf2);
  file7=strcat(dir,prefix4,suf1);
  file8=strcat(dir,prefix5,suf2);





  %%%%%%%%%%%%%%%%%%%%% READ HEADER
  %
  %
  
  numparms=3;


  fid=fopen(file1);
  [whichrnpmethod,count]=fscanf(fid,'%g',[1]);
  fclose(fid);
  

  if(whichrnpmethod==0)
    numheaders=15;
    numextras=1;
  end
  if(whichrnpmethod==1)
    numheaders=12;
    numextras=10;
  end

  
  fid=fopen(file1);
  % +1 gets whichrnpmethod again
  [myhead,count]=fscanf(fid,'%g',[numheaders+1]);
  fclose(fid);


  ii=1;
  whichrnpmethod=myhead(ii); ii=ii+1;
  ndim=myhead(ii); ii=ii+1;
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

  
  if(whichrnpmethod==0)
    % HCM and TDYN

    nhcm=myhead(ii); ii=ii+1;
    hcmmin=myhead(ii); ii=ii+1;
    hcmmax=myhead(ii); ii=ii+1;
    
    lhcmmin=log10(hcmmin);
    lhcmmax=log10(hcmmax);
    if(nhcm==1)
      steplhcm=0
    else
      steplhcm=(lhcmmax-lhcmmin)/(nhcm-1);
    end
    
    
    ntdynorye=myhead(ii); ii=ii+1;
    tdynoryemin=myhead(ii); ii=ii+1;
    tdynoryemax=myhead(ii); ii=ii+1;
    
    ltdynoryemin=log10(tdynoryemin);
    ltdynoryemax=log10(tdynoryemax);
    if(ntdynorye==1)
      stepltdynorye=0
    else
      stepltdynorye=(ltdynoryemax-ltdynoryemin)/(ntdynorye-1);
    end
    
  end

  if(whichrnpmethod==1)
    % No HCM, just rhob,T,Ye

    nhcm=1;
    hcmmin=1.0;
    hcmmax=1.0;
    
    lhcmmin=log10(hcmmin);
    lhcmmax=log10(hcmmax);
    if(nhcm==1)
      steplhcm=0
    else
      steplhcm=(lhcmmax-lhcmmin)/(nhcm-1);
    end
    
    ntdynorye=myhead(ii); ii=ii+1;
    tdynoryemin=myhead(ii); ii=ii+1;
    tdynoryemax=myhead(ii); ii=ii+1;
    
    ltdynoryemin=log10(tdynoryemin);
    ltdynoryemax=log10(tdynoryemax);
    if(ntdynorye==1)
      stepltdynorye=0
    else
      stepltdynorye=(ltdynoryemax-ltdynoryemin)/(ntdynorye-1);
    end
    
  end






  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read parameter information
  %
  %
  fid=fopen(file7);
  [myhead,count]=fscanf(fid,'%g',[numparms]);
  fclose(fid);

  ii=1;
  UTOT0=myhead(ii); ii=ii+1;
  PTOT0=myhead(ii); ii=ii+1;
  CHI0=myhead(ii); ii=ii+1;





  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % open eos.dat
  %
  %
  fid=fopen(file2);
  [mydata,count]=fscanf(fid,'%g',[nc*nrhob*ntk*nhcm*ntdynorye]);
  fclose(fid);

  % make vector into 5-D array
  temp=reshape(mydata,[nc,nrhob,ntk,nhcm,ntdynorye]);
  clear mydata;
  % set primary matricies to be each column(field) read in from SM data (i.e.
  % read data is setup with column / nx / ny / nz order
  mynewdata=permute(temp,[2,3,4,5,1]);
  clear temp;
  %
  % so now all functions are F(rhob,tk,hcm,tdynorye)
  %

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % clean-up the input data
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for p=1:nrhob
    for q=1:ntk
      for r=1:nhcm
        for s=1:ntdynorye
          nanonline=0;
          for t=1:nc
            nanonline=nanonline+isnan(mynewdata(p,q,r,s,t));
          end
          if nanonline>0
            fprintf('Found nan at p=%d q=%d r=%d s=%d\n',p,q,r,s);
            % assume NaN's occur at higher densities and just copy from lower densities
            % first 4 quantities are independent variables that shouldn't change
            for t=5:nc
              mynewdata(p,q,r,s,t) = mynewdata(p-1,q,r,s,t);
            end
            %A(isnan(A))=0
          end
        end
      end
    end
  end

  
  % {tdynorye, hcm, rhob, tk, etae, npratio, p_tot, rho_tot, s_tot, p_photon,
  % p_eleposi,p_N,p_nu}, extra1,etc.
  ii=1;
  rhob(:,:,:,:) = mynewdata(:,:,:,:,ii); ii=ii+1;
  tk = mynewdata(:,:,:,:,ii); ii=ii+1;
  if(whichrnpmethod==0)
    hcm = mynewdata(:,:,:,:,ii); ii=ii+1;
  end
  if(whichrnpmethod==1)
    % just a place holder to keep code simple
    hcm = tk*0+1;
  end
  tdynorye = mynewdata(:,:,:,:,ii); ii=ii+1;
  % below ptot,utot,stot actually don't include neutrino contribution,
  % which will be reconstructed within HARM
  ptot = mynewdata(:,:,:,:,ii); ii=ii+1;
  utot = mynewdata(:,:,:,:,ii); ii=ii+1;
  stot = mynewdata(:,:,:,:,ii); ii=ii+1;

  % rest of below don't come into detailed calculations and are just
  % other quantities to interpolate from T to UTOT
  extraii=ii;
  % instead of using extra1,extra2, etc., just directly access mynewdata
%  extra = mynewdata(:,:,:,:,extraii+extraiter); ii=ii+1;

  % COMPUTE CHI (from now on have to deal with chi separately)
  chi = utot + ptot;





  % whether to force monotonicity for functions as functions of temperature
  % otherwise EOS is non-convex and multivalued for a single internal energy
  % for the HARM code, all that really matters is that chi=u+p is single valued
  %
  if forcemonotk==1

    for p=1:nrhob
      for q=1:nhcm
        for r=1:ntdynorye
          % only need to force these to be monotonic so chi is monotonic so can have single-valued inversion
          % ptot and utot (and so chi)  need to be monotonically increasing functions of temperature in order to obtain a single temperature for a single ptot,utot,chi
          utot(p,:,q,r)=monotonize(utot(p,:,q,r));

          % entropy should be monotonic so that can be used as independent variable, but only needs to be monotonic sfor P(S) and U(S), not T(S).
          %stot(p,:,q,r)=monotonize(stot(p,:,q,r));

          if 0
            % monotonize chi and recompute ptot from monotonized chi and utot
            %chi(p,:,q,r)=monotonize(chi(p,:,q,r));
            % recompute ptot from monotonized utot and chi since it's more important to be monotonic in chi and ptot
            % leaves ptot<0 if chi chopped alot, so use old method of fixing utot and ptot
            %ptot(p,:,q,r)=chi(p,:,q,r)-utot(p,:,q,r);
          end

          if 1
            ptot(p,:,q,r)=monotonize(ptot(p,:,q,r));
            chi(p,:,q,r) = utot(p,:,q,r) + ptot(p,:,q,r);
          end
	end
      end
    end

    
  end




  if utotfix==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Temporary fix for utot since I was subtracting off rhob instead of rhob
    % c^2
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    utotnew= utot + rhob - rhob.*c.*c;

    % by itself the above generates negative utot's, so must put a lower limit
    % this wouldn't occur if did things correctly in Kaz's code

    isneg = (utotnew<=0);
    myfloor = utotnew.*0 + 1E-20;
    utot(isneg) = myfloor(isneg);

  end



  if utotdegenanalytic==1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % set 0 of utot, ptot, and chi
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %utot = utot - udegenfit_fun(rhob);
    %ptot = ptot - pdegenfit_fun(rhob);
    % chi is not read-in, and since we just modified utot and ptot, now anything (such as chi) computed from these will already be offset

    % instead of above subtraction, treat as general offset, so that output is uniform when including offset, but things are computed with standard utot and ptot
    utotoffset = udegenfit_fun(rhob);
    ptotoffset = pdegenfit_fun(rhob);
    % essentially, this means offset's only affect gridding of chosen u,p,chi as dependent variables
    % that is, we need to modify lu

    % set chidegenfit
    chidegenfit=udegenfit_fun(rhob) + pdegenfit_fun(rhob);
    chioffset  = utotoffset+ptotoffset;




  end



  if utotdegenanalytic==0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % set 0 of utot, ptot, and chi
    %		
    % Assumes monotonic utot,ptot,chi with temperature
    %		
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % pick lowest temperature as bottom
    % Duplicate this for other temperatures
    if(0)
      for q=1:ntk

        udegenfit(:,q,:,:) = utot(:,1,:,:);
        pdegenfit(:,q,:,:) = ptot(:,1,:,:);
        chidegenfit(:,q,:,:) = chi(:,1,:,:); % chi treated independently
      end
    end

    if(1)
      % don't assume monotonicity
      for p=1:ntdynorye
        for o=1:nhcm
          for m=1:nrhob

            for n=1:ntk
              mytemporary(n) = utot(m,n,o,p);
            end
            mytempmin = min(mytemporary);
            mytempmax = max(mytemporary);
            for n=1:ntk
              udegenfit(m,n,o,p) = mytempmin;
              umax(m,n,o,p) = mytempmax;
            end

            for n=1:ntk
              mytemporary(n) = ptot(m,n,o,p);
            end
            mytempmin = min(mytemporary);
            mytempmax = max(mytemporary);
            for n=1:ntk
              pdegenfit(m,n,o,p) = mytempmin;
              pmax(m,n,o,p) = mytempmax;
            end

            for n=1:ntk
              mytemporary(n) = chi(m,n,o,p);
            end
            mytempmin = min(mytemporary);
            mytempmax = max(mytemporary);
            for n=1:ntk
              % chi treated independently
              chidegenfit(m,n,o,p) = mytempmin;
              chimax(m,n,o,p) = mytempmax;
            end

          end
        end
      end
    end

    % now include multiplicative offset, where offset was chosen ahead of time so that utotdiff,ptotdiff are small and positive for all temperatures for each rhob,hcm,tdynorye
    % This procedure doesn't make sense if the degen value is near 0, then need to offset by some fraction of (max-min)
    % chi treated independently

    %utotoffset = udegenfit   - abs(udegenfit)*(1-UTOT0);
    %ptotoffset = pdegenfit   - abs(udegenfit)*(1-PTOT0);
    %chioffset  = chidegenfit - abs(udegenfit)*(1-CHI0);

    % general approach so log-intervals always equally resolved in log of temperature
    % concept is that minimum sets no scale.  Scale only set by maximum (assumes maximum is not too far from where quantity becomes as degenerate as not degenerate

    numshift = ntk;
    
    UTOTF=1.0-1E-14;
    %UTOTF=UTOT0;
    UTOTD=10.^(log10(1-UTOT0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(1-UTOTF)-log10(1-UTOT0)));
    PTOTF=1.0-1E-14;
    %PTOTF=PTOT0;
    PTOTD=10.^(log10(1-PTOT0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(1-PTOTF)-log10(1-PTOT0)));
    CHIF=1.0-1E-14;
    %CHIF=CHI0;
    CHID=10.^(log10(1-CHI0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(1-CHIF)-log10(1-CHI0)));
    
    utotoffset = udegenfit   - max(10.^(log10(abs(umax))-numshift),abs(udegenfit).*(UTOTD));
    ptotoffset = pdegenfit   - max(10.^(log10(abs(pmax))-numshift),abs(pdegenfit).*(PTOTD));
    chioffset  = chidegenfit - max(10.^(log10(abs(chimax))-numshift),abs(chidegenfit).*(CHID));


  end





  if utotdegencut==1

    utotdiff = utot - utotoffset; % used for grid of data
    ptotdiff = ptot - ptotoffset; % used for grid of data
    chidiff  = chi  - chioffset;  % used for grid of data

  end

  if utotdegencut==0

    utotdiff = utot; % used for grid of data
    ptotdiff = ptot; % used for grid of data
    chidiff  = chi;  % used for grid of data

  end







  %%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Compute sound speed from temperature-based data.  This is taken from HELM TM EOS
  % TM took it from:: references: cox & giuli chapter 24 ; timmes & swesty apj 1999
  %
  %%%%%%%%%%%%%%%%%%%%%%%

  % always set this so can compare to other sound speed
  if 1 || preinterpsoundspeed==1 ||  preinterpsoundspeed==2

    %..the temperature and density exponents (c&g 9.81 9.82) 
    %..the specific heat at constant volume (c&g 9.92)
    %..the third adiabatic exponent (c&g 9.93)
    %..the first adiabatic exponent (c&g 9.97) 
    %..the second adiabatic exponent (c&g 9.105)
    %..the specific heat at constant pressure (c&g 9.98) 
    %..and relativistic formula for the sound speed (c&g 14.29)

    %      JCM: my definitions from HELM TM EOS
    den=rhob.*c.*c;
    pres=ptot;
    temp=tk.*kb;
    % specific energy
    ener=utot./den;
    % speed of light
    %       clight=c;

    
    % get derivatives
    for q=1:nhcm
      for r=1:ntdynorye


        roughptot(:,:)=ptot(:,:,q,r);
        roughener(:,:)=ener(:,:,q,r);

        % Smooth functions since sharp boundaries cause speed to be artificially large
        myptot(:,:)=moving_average2(roughptot(:,:),1,1);
        myener(:,:)=moving_average2(roughener(:,:),1,1);
        
        % use rho c^2 so result is dimensionless
        rhocsqi=den(:,1,1,1)';
        %		rhocsqi2d=den(:,:,1,1);
        %		tsize=size(rhocsqi);
        %		sizerho=tsize(1);
        %		irho=1:sizerho;
        %		irhooffset=0.5:sizerho-0.5;
        lrhocsqi = log(rhocsqi);
        %		dlrhocsqi = gradient(lrhocsqi);

        % use kb*T so result is dimensionless
        Ti=temp(1,:,1,1)';
        %		Ti2d=temp(:,:,1,1);
        %		tsize=size(Ti);
        %		sizeT=tsize(1);
        %		iT=1:sizeT;
        %		iToffset=0.5:sizeT-0.5;
        lTi = log(Ti);
        %		dlTi = gradient(lTi);
        
        % given PofU(rho0,U,H) then derivative is dU, drho0, dH for independent variable order
        %

        % log method
        if 1
          
          [dpresdlt(:,:), dpresdld(:,:)] = gradient(myptot(:,:),lTi,lrhocsqi);
          [denerdlt(:,:), denerdld(:,:)] = gradient(myener(:,:),lTi,lrhocsqi);

          % correct for gradient position offset
          % That is, gradient takes simple difference divided by simple difference of positions
          % It does not interpolate back to central location
          % GODMARK: Apparently only at edges does it reduce to simple difference.
          % Otherwise it does put gradient at function location
          %dpresdlt(:)
          %dpresdld(:)
          %denerdlt(:)
          %denerdld(:)
          

          % Get actual derivative
          dpresdt(:,:,q,r)=dpresdlt(:,:)./temp(:,:,q,r);
          dpresdd(:,:,q,r)=dpresdld(:,:)./den(:,:,q,r);

          denerdt(:,:,q,r)=denerdlt(:,:)./temp(:,:,q,r);
          denerdd(:,:,q,r)=denerdld(:,:)./den(:,:,q,r);
        end
        % non-log method
        if 0
          [dpresdt(:,:), dpresdd(:,:)] = gradient(myptot(:,:),Ti,rhocsqi);
          [denerdt(:,:), denerdd(:,:)] = gradient(myener(:,:),Ti,rhocsqi);

        end
        

      end
    end

    
    %      dP/dT|rhob
    %       dpresdt

    %      dP/drhob|T
    %       dpresdd

    %      dener/dT|rhob
    %       denerdt

    %      HELM TM EOS extra definitions
    deni=1./den;

    %      copied from HELM TM EOS code
    zz    = pres .* deni;
    zzi   = den ./ pres;
    chit  = (temp ./ pres) .* dpresdt;
    chid  = (abs(dpresdd)+1E-20) .* zzi;
    % don't allow negative cv
    cv    = abs(denerdt)+1E-20;
    x     = (zz .* chit) ./ (temp .* cv);
    gam3  = x + 1.0d0;
    gam1  = (chit .* x) + chid;
    nabad = x ./ gam1;
    gam2  = 1.0d0 ./ (1.0d0 - nabad);
    cp    = (cv .* gam1) ./ chid;
    z     = 1.0d0 + (ener + 1.0d0) .* zzi;
    
    %      Finally quantity is sound speed squared in dimensionless units
    cs2rhoT=(gam1 ./ z);
    
    % Because of how computed, this may be negative or >1.0, so at least
    % enforce not negative
    cs2rhoT=abs(cs2rhoT)+1E-20;

  end

  

  %%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Output monotonized EOS so can check in SM
  %
  % Assumes nutotout is same size as nptotout and nchiout
  %
  % Treat chi as independent
  %
  %%%%%%%%%%%%%%%%%%%%%%%
  fid=fopen(file8,'w');
  for p=1:ntdynorye
    for o=1:nhcm
      for n=1:ntk
        for m=1:nrhob
          fprintf(fid,'%3d %3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g', ...
                  m-1, n-1, o-1, p-1, ...
                  rhob(m,n,o,p), tk(m,n,o,p), hcm(m,n,o,p), tdynorye(m,n,o,p), ...
                  ptot(m,n,o,p), utot(m,n,o,p), chi(m,n,o,p), stot(m,n,o,p),  ...
                  pdegenfit(m,n,o,p), udegenfit(m,n,o,p), chidegenfit(m,n,o,p), ...
                  ptotoffset(m,n,o,p), utotoffset(m,n,o,p), chioffset(m,n,o,p), ...
                  ptotdiff(m,n,o,p), utotdiff(m,n,o,p), chidiff(m,n,o,p), ...
                  cs2rhoT(m,n,o,p) ...
                  );
          for ei=extraii:extraii+numextras-1
            fprintf(fid,'%21.15g ', mynewdata(m,n,o,p,ei));
          end
          %                fprintf(fid,'%3d %3d %3d %3d', ...
          %                m-1, n-1, o-1, p-1 ...
          %		);
          
          fprintf(fid,'\n');
        end
      end
    end
  end
  fclose(fid);







  % DEBUG
  % for ii=1:48 fprintf('%21.15g %21.15g %21.15g\n',log10(utot(28,ii,1,1)/(c*c)),log10(utotnew(28,ii,1,1)/(c*c)),log10(tk(28,ii,1,1))); end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % some computed things
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  rhobcsq		= rhob.*c.*c;
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
    [mydata,count]=fscanf(fid,'%g',[nco*nrhob*ntk*nhcm*ntdynorye]);
    fclose(fid);

    % make vector into 5-D array
    temp=reshape(mydata,[nco,ntk,nrhob,nhcm,ntdynorye]);
    clear mydata;
    % set primary matricies to be each column(field) read in from SM data (i.e.
    % read data is setup with column / nx / ny / nz order
    mynewdata=permute(temp,[2,3,4,5,1]);
    clear temp;
    %
    %
    % so now all functions are F(rhob,tk,hcm,tdynorye)
    %

    %


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
  % MARKER:
  % Using temperature as marker to tell if within original Kaz EOS or outside and within
  % interpolated/extrapolted/filled-in region
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % First construct 1-D grid of new dependent quantities
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % factor by which to enhance resolution of internal values before downgrading to output resolution
  factor=10;
  %factor=1;



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % grid of internal energy (utot) (used to be also grid for pressure (ptot), and \chi=u+p)
  nutot=ntk*factor;
  %lutotmin=log10(rhobmin*c*c);
  %lutotmax=log10(rhobmax*c*c);

  if utotfix==0
    lutotmin=log10(min(min(min(min(utotdiff)))));
    lutotmax=log10(max(max(max(max(utotdiff)))));
  end
  if utotfix==1
    % presently utot is wrong after fixup when small due to rhob subraction
    % choose reasonable lower limit
    lutotmin=15.0;
    lutotmax=log10(max(max(max(max(utotdiff)))));
  end

  steplutot = (lutotmax-lutotmin)/(nutot-1);
  lutotgrid=lutotmin:steplutot:lutotmax;
  utotgrid = 10.^lutotgrid;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % grid of ptot
  nptot=nutot;
  lptotmin=log10(min(min(min(min(ptotdiff)))));
  lptotmax=log10(max(max(max(max(ptotdiff)))));
  steplptot = (lptotmax-lptotmin)/(nptot-1);
  lptotgrid=lptotmin:steplptot:lptotmax;
  ptotgrid = 10.^lptotgrid;



  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % grid of \chi=u+p
  nchi=nutot;
  lchimin=log10(min(min(min(min(chidiff)))));
  lchimax=log10(max(max(max(max(chidiff)))));
  steplchi = (lchimax-lchimin)/(nchi-1);
  lchigrid=lchimin:steplchi:lchimax;
  chigrid = 10.^lchigrid;



  % entropy as a variable is only used temporarily
  nsspec=nutot;
  lsspecmin=log10(min(min(min(min(sspec)))));
  lsspecmax=log10(max(max(max(max(sspec)))));
  steplsspec = (lsspecmax-lsspecmin)/(nsspec-1);
  lsspecgrid=lsspecmin:steplsspec:lsspecmax;
  sspecgrid = 10.^lsspecgrid;


  % OUTPUT grid of functions of utot
  nutotout=ntk; % anything can be put here, but assume basically same as size of ntk
                %lutotmin=log10(rhobmin*c*c);
                %lutotmax=log10(rhobmax*c*c);
  lutotoutmin=lutotmin;
  lutotoutmax=lutotmax;
  steplutotout = (lutotoutmax-lutotoutmin)/(nutotout-1);
  lutotoutgrid=lutotoutmin:steplutotout:lutotoutmax;
  utotoutgrid = 10.^lutotoutgrid;

  % OUTPUT grid of functions of ptot
  nptotout=nutotout;
  %lptotmin=log10(rhobmin*c*c);
  %lptotmax=log10(rhobmax*c*c);
  lptotoutmin=lptotmin;
  lptotoutmax=lptotmax;
  steplptotout = (lptotoutmax-lptotoutmin)/(nptotout-1);
  lptotoutgrid=lptotoutmin:steplptotout:lptotoutmax;
  ptotoutgrid = 10.^lptotoutgrid;

  % OUTPUT grid of functions of chi
  nchiout=nutotout;
  %lchimin=log10(rhobmin*c*c);
  %lchimax=log10(rhobmax*c*c);
  lchioutmin=lchimin;
  lchioutmax=lchimax;
  steplchiout = (lchioutmax-lchioutmin)/(nchiout-1);
  lchioutgrid=lchioutmin:steplchiout:lchioutmax;
  chioutgrid = 10.^lchioutgrid;





  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Now interpolate to obtain desired functions
  %
  % interpolate in log10
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  clear extraofU;

  
  % say X, Y, Xi and get Yi
  for p=1:nrhob
    for q=1:nhcm
      for r=1:ntdynorye

        % p
        %q
        %r

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % BEGIN OLD CODE
        %
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

        %
        % END OLD CODE
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        % make sure things interpolating are unique
        % Here x = utotdiff,ptotdiff,chidiff,sspec and y = ptot, hspec,
        % utot, ptot, stot, Tk, extra?, etc.
        

        [lutotptotx,lutotptoty] = consolidator(log10(utotdiff(p,:,q,r)),log10(ptot(p,:,q,r)),'mean',CONTOL);

        [lhspecx,lhspecy] = consolidator(log10(utotdiff(p,:,q,r)),log10(hspec(p,:,q,r)),'mean',CONTOL);

        [lcs2rhoTx,lcs2rhoTy] = consolidator(log10(utotdiff(p,:,q,r)),log10(cs2rhoT(p,:,q,r)),'mean',CONTOL);

        [lptotutotx,lptotutoty] = consolidator(log10(ptotdiff(p,:,q,r)),log10(utot(p,:,q,r)),'mean',CONTOL);

        [lchiptotx,lchiptoty] = consolidator(log10(chidiff(p,:,q,r)),log10(ptot(p,:,q,r)),'mean',CONTOL);

        [lsspecptotx,lsspecptoty] = consolidator(log10(sspec(p,:,q,r)),log10(ptot(p,:,q,r)),'mean',CONTOL);
        [lsspecutotx,lsspecutoty] = consolidator(log10(sspec(p,:,q,r)),log10(utot(p,:,q,r)),'mean',CONTOL);

        [lutotstotx,lutotstoty] = consolidator(log10(utotdiff(p,:,q,r)),log10(stot(p,:,q,r)),'mean',CONTOL);

        % Temperature is used to establish if really within valid EOS (there was a valid inversion)

        % Tk(u)
        [lutottkx,lutottky] = consolidator(log10(utotdiff(p,:,q,r)),log10(tk(p,:,q,r)),'mean',CONTOL);

        % Tk(p)
        [lptottkx,lptottky] = consolidator(log10(ptotdiff(p,:,q,r)),log10(tk(p,:,q,r)),'mean',CONTOL);

        % Tk(chi)
        [lchitkx,lchitky] = consolidator(log10(chidiff(p,:,q,r)),log10(tk(p,:,q,r)),'mean',CONTOL);


        % EXTRAS
        % mynewdata(:,:,:,:,ei)
        % assume extras are all >0 (replace 0 with 1E-20)
        %p
        %  q
        %  r
        %  size(lutotextrax)
        %  size(utotdiff)
        %  size(mynewdata)
        
        %setup size of extras array
        %lutotextrax = zeros(numextras,ntk);
        %lutotextray = zeros(numextras,ntk);


        [tempx tempy] = consolidator(log10(utotdiff(p,:,q,r)),log10(abs(mynewdata(p,:,q,r,ei))+1E-20),'mean',CONTOL);
        sizex=size(tempx(:,1));
        sizexx=sizex(1);
        lutotextrax=zeros(numextras,sizexx);
        sizey=size(tempy(:,1));
        sizeyy=sizey(1);
        lutotextray=zeros(numextras,sizeyy);

        for ei=extraii:extraii+numextras-1
          % for extra : comes second for stupid reason that
          % size(utotdiff(p,:,q,r)) is 1 200
          [lutotextrax(ei-extraii+1,:) lutotextray(ei-extraii+1,:)] = consolidator(log10(utotdiff(p,:,q,r)),log10(abs(mynewdata(p,:,q,r,ei))+1E-20),'mean',CONTOL);
        end
        % now transpose result so that 200 1
        temp=lutotextrax';
        lutotextrax=temp;
        temp=lutotextray';
        lutotextray=temp;

        

        %%%%%%% BEGIN DEBUG
        %if		p==197
        
        %    lsspecptotx'
        %    lsspecutoty'


        % for mmm=1:200 fprintf('%21.15g\n',sspec(p,mmm,q,r)); end
        %lutotptotx'
        %lutotptoty'

        %end

        %  f(x) = interp(x(t),f(t),x_i)
        %%%%%%% END DEBUG
        


        %%%%%%%%%%%%%%%% Interpolate F(rho0,?)
        %
        % BELOW dependent variables are actually differenced (offsetted) versions
        % below is P(rho0,u,H)
        % below directly used as P(rho0,u)
        % dP/du |rho0
        % dP/drho0 | u
        PofU(p,:,q,r) = 10.^(myinterp1(1,lutotptotx, lutotptoty, lutotgrid','linear'));

        HofU(p,:,q,r) = 10.^(myinterp1(2,lhspecx, lhspecy, lutotgrid','linear'));

        % c_s^2 in dimensionless units as function of rhob,U
        cs2(p,:,q,r) = 10.^(myinterp1(2,lcs2rhoTx, lcs2rhoTy, lutotgrid','linear'));

        % below is u(rho0,P,H)
        % below directly used as u(rho0,p)
        UofP(p,:,q,r) = 10.^(myinterp1(3,lptotutotx, lptotutoty, lptotgrid','linear'));

        % below is P(rho0,\chi,H)
        % below used for P[rho0,chi]
        % 1/ (dchi/dp)|rho0
        % 1/(drho0/dp)|chi
        PofCHI(p,:,q,r) = 10.^(myinterp1(4,lchiptotx, lchiptoty, lchigrid','linear'));

        % Below is PofS(rho0,S,H)
        % below used for c_s^2 = 1/h dp/drho0|S
        PofS(p,:,q,r) = 10.^(myinterp1(5,lsspecptotx, lsspecptoty, lsspecgrid','linear'));
        UofS(p,:,q,r) = 10.^(myinterp1(6,lsspecptotx, lsspecutoty, lsspecgrid','linear'));

        % Below is SofU(rho0,u,H)
        SofU(p,:,q,r) = 10.^(myinterp1(8,lutotstotx, lutotstoty, lutotgrid','linear'));

        % Below is T(rho0,u)
        tkofU(p,:,q,r) = 10.^(myinterp1(29,lutottkx, lutottky, lutotgrid','linear'));

        % Below is T(rho0,p)
        tkofP(p,:,q,r) = 10.^(myinterp1(30,lptottkx, lptottky, lptotgrid','linear'));

        % Below is T(rho0,chi)
        tkofCHI(p,:,q,r) = 10.^(myinterp1(31,lchitkx, lchitky, lchigrid','linear'));

        % Below is extra1ofU(rho0,u,H)
        % mynewdata(:,:,:,:,ei)
        % so far all these things are positive definite, so log interpolate ok
        for ei=1:numextras
          % ei comes first for issues of size()
          extraofU(p,:,q,r,ei) = 10.^(myinterp1(7,lutotextrax(:,ei), lutotextray(:,ei), lutotgrid','linear'));
        end

        
        
      end
    end
  end




  %DEBUG
  %HofU(1,14,1,1)
  %HofU(1,:,1,1)
  %DEBUG






  if doclean
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Clean/fixup F(rho0,?)
    %
    % Adjust quantities to be physical when interpolated range is beyond
    % existing range
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % rhobgrid, utotgrid, ptotgrid, chigrid, and sspecgrid need to be turned into full 4-D quantity for cleanvar

    % say X, Y, Xi and get Yi

    % these are super-sampled versions so far
    for p=1:nutot
      for q=1:nhcm
        for r=1:ntdynorye
          rhobgrid4d(:,p,q,r)    = rhob(:,1,1,1);
          rhobcsqgrid4d(:,p,q,r) = rhobcsq(:,1,1,1);
        end
      end
    end
    for p=1:nrhob
      for q=1:nhcm
        for r=1:ntdynorye
          utotgrid4d(p,:,q,r)  = utotgrid(:);
          ptotgrid4d(p,:,q,r)  = ptotgrid(:);
          chigrid4d(p,:,q,r)   = chigrid(:);
          sspecgrid4d(p,:,q,r) = sspecgrid(:);
        end
      end
    end

    %%%%%%%%%%% Clean F(rho0,?)
    

    PofU = cleanvar(1, PofU, rhobgrid4d, utotgrid4d);
    HofU = cleanvar(2, HofU, rhobcsqgrid4d, utotgrid4d); % specific enthalpy is dimensionless
    UofP = cleanvar(3, UofP, rhobgrid4d, ptotgrid4d);
    PofCHI = cleanvar(4, PofCHI, rhobgrid4d, chigrid4d);
    PofS = cleanvar(5, PofS, rhobcsqgrid4d, sspecgrid4d); % here S is specific entropy and uses rhobcsq
    UofS = cleanvar(6, UofS, rhobcsqgrid4d, sspecgrid4d); % here S is specific entropy and uses rhobcsq
    SofU = cleanvar(8, SofU, rhobcsqgrid4d, utotgrid4d); % here S is entropy density and uses rhobcsq

    %size(tkofU)
    %size(rhobgrid4d)
    %size(utotgrid4d)

    tkofU = cleanvar(29, tkofU, rhobgrid4d, utotgrid4d);
    tkofP = cleanvar(30, tkofP, rhobgrid4d, ptotgrid4d);
    tkofCHI = cleanvar(31, tkofCHI, rhobgrid4d, chigrid4d);

    cs2 = cleanvar(40, cs2, rhobcsqgrid4d, utotgrid4d); % here cs2 is dimensionless

    % extras:
    % mynewdata(:,:,:,:,ei)
%    for ei=1:numextras
      %extraofU(ei,:,:,:,:) = cleanvar(7, extraofU(ei,:,:,:,:),rhobgrid4d, utotgrid4d);
    % clean all at once
    extraofU = cleanvar(7, extraofU, rhobgrid4d, utotgrid4d);
%    end



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
  % stoti = utotgrid;
  % dstoti=gradient(stoti);

  % specific entropy as independent quantity
  % sspeci = utotgrid./rhocsqi;
  sspeci = sspecgrid;
  dsspeci=gradient(sspeci);

  utoti = utotgrid;
  dutoti=gradient(utoti);

  ptoti = ptotgrid;
  dptoti=gradient(ptoti);

  chii = chigrid;
  dchii=gradient(chii);


  
  % say X, Y, Xi and get Yi
  for q=1:nhcm
    for r=1:ntdynorye
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





  %%%%%%%%%%%%%%%%%%
  %
  % Done with temporary variables used for gradient operator
  %
  %%%%%%%%%%%%%%%%%%

  clear rhocsqi drhocsqi rhoi drhoi sspeci dsspeci utoti dutoti ptoti dptoti chii dchii




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Some derivatives need to be have a change of variable
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






  % say X, Y, Xi and get Yi
  for p=1:nrhob
    for q=1:nhcm
      for r=1:ntdynorye
        
        % below is presently dP/drho(rho0,S)|S and need dP/drho(rho0,u)|S so change S->u
        % dPofSdrho0(:,:,q,r)

        % GODMARK: Should below be utot(rho0,T) -> UofS(rho0,S) ?  Seems so.

        [lutotdpx,lutotdpy] = consolidator(log10(UofS(p,:,q,r)),dPofSdrho0(p,:,q,r),'mean',CONTOL);

        % Below is dP/drho0(U,rho0,H)|S
        % below used for c_s^2 = 1/h dp/drho0|S
        dPofUdrho0cS(p,:,q,r) = myinterp1(9,lutotdpx, lutotdpy, lutotgrid','linear');

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

    dPofUdrho0cS = cleanvar(9, dPofUdrho0cS, rhobcsqgrid4d, utotgrid4d);

  end




  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Now compute things that don't involve derivatives
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if preinterpsoundspeed==0 || preinterpsoundspeed==2 

    % c_s^2 in (cm/s)^2
    % dPofUdrho0cS already has 1/c^2 in it
    % so below is per unit c^2
    cs2post = (1.0./HofU).*(dPofUdrho0cS);

  end

  if preinterpsoundspeed==0

    cs2=cs2post;

  end

  if preinterpsoundspeed==2
    % Use cs2post to fix cs2
    for p=1:ntdynorye
      for o=1:nhcm
        for n=1:nutot % supersampled so far
          for m=1:nrhob
            
            if cs2(m,n,o,p)>2.0*cs2post(m,n,o,p) && cs2post(m,n,o,p)>=0.0
              cs2(m,n,o,p) = cs2post(m,n,o,p);
            end
            if cs2(m,n,o,p)<=0.0 && cs2post(m,n,o,p)>0.0
              if cs2post(m,n,o,p)<1.0
                cs2(m,n,o,p) = cs2post(m,n,o,p);
              else
                cs2(m,n,o,p) = 1E-30;
              end
            end

          end
        end
      end
    end



  end


  % do in any case
  cs2cgs = cs2.*c.*c;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %		list of things to output
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


  % extraofU(rho0,u,H,ei)

  % tkofU(rho0,u,H)




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
      for r=1:ntdynorye
        
        %            [lutotdpx,lutotdpy] = consolidator(log10(UofS(p,:,q,r)),log10(dPofSdrho0(p,:,q,r)),'mean',CONTOL);

        PofUout(p,:,q,r)          = 10.^(myinterp1(1,lutotgrid, log10(PofU(p,:,q,r)), lutotoutgrid','linear'));
        HofUout(p,:,q,r)          = 10.^(myinterp1(2,lutotgrid, log10(HofU(p,:,q,r)), lutotoutgrid','linear'));
        UofPout(p,:,q,r)          = 10.^(myinterp1(3,lptotgrid, log10(UofP(p,:,q,r)), lptotoutgrid','linear'));
        % dPofUdrho0out(p,:,q,r)   = 10.^(myinterp1(20,lutotgrid, log10(dPofUdrho0(p,:,q,r)), lutotoutgrid','linear'));
        % dPofUduout(p,:,q,r)      = 10.^(myinterp1(21,lutotgrid, log10(dPofUdu(p,:,q,r)), lutotoutgrid','linear'));
        dPofUdrho0out(p,:,q,r)    = myinterp1(20,lutotgrid, dPofUdrho0(p,:,q,r), lutotoutgrid','linear');
        dPofUduout(p,:,q,r)       = myinterp1(21,lutotgrid, dPofUdu(p,:,q,r), lutotoutgrid','linear');
        cs2out(p,:,q,r)           = myinterp1(22,lutotgrid, cs2(p,:,q,r), lutotoutgrid','cubic');
        cs2cgsout(p,:,q,r)        = myinterp1(23,lutotgrid, cs2cgs(p,:,q,r), lutotoutgrid','cubic');
        SofUout(p,:,q,r)          = 10.^(myinterp1(8,lutotgrid, log10(SofU(p,:,q,r)), lutotoutgrid','linear'));
        % dSofUdrho0out(p,:,q,r)   = 10.^(myinterp1(24,lutotgrid, log10(dSofUdrho0(p,:,q,r)), lutotoutgrid','linear'));
        % dSofUduout(p,:,q,r)      = 10.^(myinterp1(25,lutotgrid, log10(dSofUdu(p,:,q,r)), lutotoutgrid','linear'));
        dSofUdrho0out(p,:,q,r)    = myinterp1(24,lutotgrid, dSofUdrho0(p,:,q,r), lutotoutgrid','linear');
        dSofUduout(p,:,q,r)       = myinterp1(25,lutotgrid, dSofUdu(p,:,q,r), lutotoutgrid','linear');
        PofCHIout(p,:,q,r)        = 10.^(myinterp1(4,lchigrid, log10(PofCHI(p,:,q,r)), lchioutgrid','linear'));
        % dPofCHIdrho0out(p,:,q,r) = 10.^(myinterp1(26,lchigrid, log10(dPofCHIdrho0(p,:,q,r)), lchioutgrid','linear'));
        % dPofCHIdchiout(p,:,q,r)  = 10.^(myinterp1(27,lchigrid, log10(dPofCHIdchi(p,:,q,r)), lchioutgrid','linear'));
        dPofCHIdrho0out(p,:,q,r)  = myinterp1(26,lchigrid, dPofCHIdrho0(p,:,q,r), lchioutgrid','linear');
        dPofCHIdchiout(p,:,q,r)   = myinterp1(27,lchigrid, dPofCHIdchi(p,:,q,r), lchioutgrid','linear');
        
        tkofUout(p,:,q,r)         = 10.^(myinterp1(29,lutotgrid, log10(tkofU(p,:,q,r)), lutotoutgrid','linear'));
        tkofPout(p,:,q,r)         = 10.^(myinterp1(30,lptotgrid, log10(tkofP(p,:,q,r)), lptotoutgrid','linear'));
        tkofCHIout(p,:,q,r)       = 10.^(myinterp1(31,lchigrid, log10(tkofCHI(p,:,q,r)), lchioutgrid','linear'));

        % extras:
        % mynewdata(:,:,:,:,ei)
        for ei=1:numextras
          extraofUout(p,:,q,r,ei) = 10.^(myinterp1(28,lutotgrid, log10(extraofU(p,:,q,r,ei)), lutotoutgrid','linear'));
        end
        
        
      end
    end
  end

  %HofUout(1,:,1,1)



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Clear interpolated but pre-output version of variables
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  clear HofU cs2 PofU UofP dPofUdrho0 dPofUdu cs2cgs SofU dSofUdrho0 dSofUdu PofCHI dPofCHIdrho0 dPofCHIdchi tkofU tkofP tkofCHI;
  clear extraofU;



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Adjust quantities to be physical in case of numerical error
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  fprintf('Min c_s^2/c^2\n');
  [mincs2out,mincs2outI]=min(min(min(min(cs2out))))
  fprintf('Max c_s^2/c^2\n');
  [maxcs2out,maxcs2outI]=max(max(max(max(cs2out))))

  fprintf('Min c_s^2[cgs]\n');
  [mincs2cgsout,mincs2cgsoutI]=min(min(min(min(cs2cgsout))))
  fprintf('Max c_s^2[cgs]\n');
  [maxcs2cgsout,maxcs2cgsoutI]=max(max(max(max(cs2cgsout))))

  % adjust speed of sound to be no smaller than 0 and no larger than 1

  for p=1:ntdynorye
    for o=1:nhcm
      for n=1:nutotout
        for m=1:nrhob
          
          if cs2out(m,n,o,p)>1.0
            cs2out(m,n,o,p)=1.0-CONTOL;
          end
          if cs2out(m,n,o,p)<0.0
            cs2out(m,n,o,p)=0.0;
          end

          if cs2cgsout(m,n,o,p)>c.*c
            cs2cgsout(m,n,o,p)=(1.0-CONTOL).*c.*c;
          end
          if cs2cgsout(m,n,o,p)<0.0
            cs2cgsout(m,n,o,p)=0.0;
          end

        end
      end
    end
  end

  fprintf('New Min c_s^2/c^2\n');
  [mincs2out,mincs2outI]=min(min(min(min(cs2out))))
  fprintf('New Max c_s^2/c^2\n');
  [maxcs2out,maxcs2outI]=max(max(max(max(cs2out))))

  fprintf('New Min c_s^2[cgs]\n');
  [mincs2cgsout,mincs2cgsoutI]=min(min(min(min(cs2cgsout))))
  fprintf('New Max c_s^2[cgs]\n');
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
  % Notice that below indicies are in C-order and C-style (start with 0 instead of 1)
  %
  % Assumes nutotout is same size as nptotout and nchiout
  %
  %%%%%%%%%%%%%%%%%%%%%%%

  %number of columns outputted (including numbers indicating position of element in table)
  %NUMOUTCOLUMNS=25 % was 23 before new degen offset method
  NUMOUTCOLUMNS=4+1+3+2+2+2+1+3+3+3+numextras;


  fid=fopen(file3,'w');
  for p=1:ntdynorye
    for o=1:nhcm
      for n=1:nutotout
        for m=1:nrhob
          fprintf(fid,'%3d %3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g', ...
                  m-1, n-1, o-1, p-1, ...
                  rhob(m,n,o,p), utotoutgrid(n), ptotoutgrid(n), chioutgrid(n), hcm(m,n,o,p), tdynorye(m,n,o,p), ...
                  PofUout(m,n,o,p), UofPout(m,n,o,p), ...
                  dPofUdrho0out(m,n,o,p), dPofUduout(m,n,o,p), ...
                  cs2cgsout(m,n,o,p), ...
                  SofUout(m,n,o,p), dSofUdrho0out(m,n,o,p), dSofUduout(m,n,o,p), ...
                  PofCHIout(m,n,o,p), dPofCHIdrho0out(m,n,o,p), dPofCHIdchiout(m,n,o,p), ...
                  tkofUout(m,n,o,p),tkofPout(m,n,o,p),tkofCHIout(m,n,o,p) ...
                  );
          % mynewdata(:,:,:,:,ei)
          for ei=1:numextras
            fprintf(fid,'%21.15g ',extraofUout(m,n,o,p,ei));
          end
          fprintf(fid,'\n');
        end
      end
    end
  end
  fclose(fid);







  %%%%%%%%%%%%%%%%%%%%%%%%
  %
  %
  %  Output utot,ptot,chi as functions of rhob for T=0 (hence TDYNORYE and HCM not important -- so becomes 1D)
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%
  fid=fopen(file6,'w');

  % corresponds to n=1 solution since u~0 implies T~0.  Reduced to degenerate solution independent of temperature
  n=1;
  for p=1:ntdynorye
    for o=1:nhcm
      for m=1:nrhob
        fprintf(fid,'%3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g ', ...
		m-1, o-1, p-1, ...
		rhob(m,n,o,p), hcm(m,n,o,p), tdynorye(m,n,o,p), ...
		utotoffset(m,n,o,p), ptotoffset(m,n,o,p), chioffset(m,n,o,p) ...
		);
        fprintf(fid,'\n');
      end
    end
  end
  fclose(fid);







  % some debug stuff:
  %
  % for ii=1:48 fprintf('%21.15g %21.15g\n',log10(utot(28,ii,1,1)),log10(tk(28,ii,1,1)));end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % write header to file
  %
  % Note that steps and all such things are computed same as in Kaz's code
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % first is table size, then as read-in by HARM code:
  % second [4] : 0 = lower log_base limit, 1 = upper log_base limit, 2=step, 3 = divisor of grid position 4=base of log, 5 = linear value of offset for log_base stepping so can control how resolved


  % set base and linear offset here since not using this in eos_extract.m yet
  baselrhob=10.0;
  linearoffsetlrhob=0.0;

  baselutotout=10.0;
  linearoffsetlutotout=0.0;
  baselptotout=10.0;
  linearoffsetlptotout=0.0;
  baselchiout=10.0;
  linearoffsetlchiout=0.0;

  baselhcm=10.0;
  linearoffsetlhcm=0.0;
  baseltdynorye=10.0;
  linearoffsetltdynorye=0.0;
  baseltk=10.0;
  linearoffsetltk=0.0;

  fid=fopen(file5,'w');
  fprintf(fid,'%d %d\n',whichrnpmethod,NUMOUTCOLUMNS);	% number of output colums
  fprintf(fid,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',nrhob,lrhobmin,lrhobmax,steplrhob,baselrhob,linearoffsetlrhob);

  fprintf(fid,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',nutotout,lutotoutmin,lutotoutmax,steplutotout,baselutotout,linearoffsetlutotout);
  fprintf(fid,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',nptotout,lptotoutmin,lptotoutmax,steplptotout,baselptotout,linearoffsetlptotout);
  fprintf(fid,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',nchiout,lchioutmin,lchioutmax,steplchiout,baselchiout,linearoffsetlchiout);

  fprintf(fid,'%d	%21.15g %21.15g %21.15g %21.15g %21.15g\n',nhcm,lhcmmin,lhcmmax,steplhcm,baselhcm,linearoffsetlhcm);
  fprintf(fid,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',ntdynorye,ltdynoryemin,ltdynoryemax,stepltdynorye,baseltdynorye,linearoffsetltdynorye);
  fprintf(fid,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',ntk,ltkmin,ltkmax,stepltk,baseltk,linearoffsetltk);
  fclose(fid);
















  fclose('all');

  % BELOW FOR relativity
  if onrelativity
    quit;
  end


end






















% See phivsr.m
% should not be function of temperature (tk)
function out = udegenfitOLD(rho0)

K = 1.22929981645031E15;
A = 1.51663;
B = 4.6E6;
C = 0.9999;
D = 1.0;
Ye = yefit(rho0,tk);

out = A.*K.*((rho0+B).*Ye).^(C/3).*rho0.^D;

% from SM:
% set udegenfit1=1.51663*K*((rhob+4.6E6)*yefit)**(0.9999/3)*rhob**(1.0)

end


% See phivsr.m
% should not be function of temperature (tk)
function out = udegenfit_fun(rhob)

Ye = yefit(rhob);
K = 1.22929981645031E15;

A = 1.516;
B = 0.0E6;
C = 4.0/3.0;
D = 1.0/3.0;
udegenfit1 = A.*K.*(rhob+B).^(C).*Ye.^(D);

A = 1.66*1E2*1.516;
B = 0.0E6;
C = 1.0;
D = 1.0/3.0;
udegenfit2 = A.*K.*(rhob+B).^(C).*Ye.^(D);

npow=1.95;

% larger function wins
out = (udegenfit1.^npow + udegenfit2.^npow).^(1/npow);


%		FROM SM:
% try log with break
%		set npow=1.95
%		set udegenfit1=1.516*K*((rhob+0.0E6))**(4.0/3)* yefit**(1.0/3)
%		set udegenfit2=1.66*1E2*1.516*K*((rhob+0.0E6))**(1.0)* yefit**(1.0/3)
%		#set udegenfit = (1.0/(1.0/udegenfit1**npow+1.0/udegenfit2**npow))**(1/npow)
%		set udegenfit = (udegenfit1**npow + udegenfit2**npow)**(1/npow)

end








% See phivsr.m
function out = pdegenfit_fun(rhob)

K = 1.22929981645031E15;
A1 = 3.0310;
A2 = 0.0235;
B = 0.0;
Ye = yefit(rhob);

pdegenfit1=(A1.*K.*((rhob+B).*yefit).^(1/3).*rhob).*(4.0/3.0-1.0).*yefit;
pdegenfit2=(A2.*K.*((rhob+B).*yefit).^(2/3).*rhob).*(4.0/3.0-1.0).*yefit;
npow=1.94;

out = (1.0./(1.0./pdegenfit1.^npow+1.0./pdegenfit2.^npow)).^(1./npow);

% from SM:
% below best fit so far for ptot
%	set pdegenfit1=(3.0310*K*((rhob+0E6)*yefit)**(1/3)*rhob)*(4.0/3.0-1.0)*yefit
%		set pdegenfit2=(0.0235*K*((rhob+0E6)*yefit)**(2.00/3)*rhob)*(4.0/3.0-1.0)*yefit
%		set npow=1.94
%		set pdegenfit = (1.0/(1.0/pdegenfit1**npow+1.0/pdegenfit2**npow))**(1/npow)


end


% \chi\equiv u+p
function out = chidegenfit_fun(rhob)

out = udegenfit_fun(rhob) + pdegenfit_fun(rhob);

end






% From jon_helm.f
%     From SM, compute Ye for presupernova-like abundances
function out = yefit(rhob)

tk=0.0;

if(tk<10.^(9.532))
  yefit = 0.5;
else
  if(tk>10.^(9.92))
    yefit = 0.43;
  else
    yefit = (0.43-0.5)/(10^(9.92)-10^(9.532))*(tk-10.^(9.532))+0.5;
  end
end

out = yefit;

end


