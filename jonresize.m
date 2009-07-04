function Ilow=jonresize(Ihigh,newsize);

nylow=newsize(2)
nxlow=newsize(1)
[nyhigh nxhigh] = size(Ihigh)

Ilow=zeros(nylow,nxlow);
IlowW=zeros(nylow,nxlow);

  for jh=1:nyhigh
    for ih=1:nxhigh
      ilfrac=ih*nxlow/nxhigh;
      jlfrac=jh*nylow/nyhigh;
      ioldbcl=(ceil(ilfrac)-1)*nxhigh/nxlow;
      joldbcl=(ceil(jlfrac)-1)*nyhigh/nylow;
      deltaih=abs(ioldbcl-(ih-0.5));
      deltajh=abs(joldbcl-(jh-0.5));
      if(deltaih>1.0) 
            deltaih=1.0; 
            end
      if(deltajh>1.0) 
            deltajh=1.0; 
            end
      
      W=deltaih*deltajh;
      
      il=ceil(ilfrac);
      jl=ceil(jlfrac);
      
      Ilow(jl,il)=Ilow(jl,il)+Ihigh(jh,ih)*W;
      IlowW(jl,il)=IlowW(jl,il)+W;
      
      %if((il==1)&(jl==1))
      %ilfrac
      %jlfrac
      %ioldbcl
      %joldbcl
      %deltaih
      %deltajh
      %W
      %il
      %jl
      %Ihigh(jh,ih)
      %Ilow(jl,il)
      %IlowW(jl,il)
      end
    end    
    end
  
  for jl=1:nylow
    for il=1:nxlow
      Ilow(jl,il)=Ilow(jl,il)/IlowW(jl,il);
  end
 end


