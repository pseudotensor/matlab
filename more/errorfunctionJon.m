function errorfunctionJon(i0, i1, m)
mminmax =[

    0.8560    1.1069
   -0.3420    0.3420
   -6.1667    7.9483
   -0.3432    0.3432
    0.8560    1.1105
   -7.4616    8.7055
    0.8500    1.1637
   -5.9350    6.8660];


%plots error function for one paramters (others being held at global min)
%mrest: 7 fized parameters
[r c s]=size(i0); %r=ny; c=nx; s=nz
[rp cp sp]=size(i1);

% jon's 1/f change
%i0(find(abs(i0)<.0001))=sign(i0(find(abs(i0)<.0001)))*.0001;
%i1(find(abs(i1)<.0001))=sign(i1(find(abs(i1)<.0001)))*.0001;
%i0=1/i0
%i1=1/i1
% end jon's 1/f change

i1vect=zeros(rp*cp*sp, 4);

%vectorize i0
i0vect=zeros(r*c*s, 4); % vectorized i0 has 4 columns - with y, x, z coordinates and intensity value in the fourth column)
h=1:r*c*s;
i0vect(:, 1)=(mod(h-1, r)+1)';%y
i0vect(:, 2)=(floor(mod((h-1)/r, c))+1)'; %x
i0vect(:, 3)=(floor((h-1)/(c*r))+1)'; %z
i0vect(:, 4) = i0(:);



numi=100
listpari=1:numi
for mytry = 1:8
for par = 0:numi
    if(mytry==1)
        di=4/numi
        si=-2
        pari=si+di*par
        mmytry=[pari m(2) m(3) m(4) m(5) m(6) m(7) m(8)];
    end
    if(mytry==2)
        di=4/numi
        si=-2
        pari=si+di*par
        mmytry=[m(1) pari m(3) m(4) m(5) m(6) m(7) m(8)];
    end
    if(mytry==3)
        di=2*r/numi
        si=-r
        pari=si+di*par
        mmytry=[m(1) m(2) pari m(4) m(5) m(6) m(7) m(8)];
    end
    if(mytry==4)
        di=4/numi
        si=-2
        pari=si+di*par
        mmytry=[m(1) m(2) m(3) pari m(5) m(6) m(7) m(8)];
    end
    if(mytry==5)
        di=4/numi
        si=-2
        pari=si+di*par
        mmytry=[m(1) m(2) m(3) m(4) pari m(6) m(7) m(8)];
    end
    if(mytry==6)
        di=2*c/numi
        si=-c
        pari=si+di*par
        mmytry=[m(1) m(2) m(3) m(4) m(5) pari m(7) m(8)];
    end
    if(mytry==7)
        di=4/numi
        si=-2
        pari=si+di*par
        mmytry=[m(1) m(2) m(3) m(4) m(5) m(6) pari m(8)];
    end
    if(mytry==8)
        di=2*s/numi
        si=-s
        pari=si+di*par
        mmytry=[m(1) m(2) m(3) m(4) m(5) m(6) m(7) pari];
    end
% i1vect(:, 1)=m(4)*i0vect(:, 2)+m(5)*i0vect(:, 1)+m(6); %yp
%    i1vect(:, 2)=m(1)*i0vect(:, 2)+m(2)*i0vect(:, 1)+m(3);%xp
%    i1vect(:, 3)=m(7)*i0vect(:, 3) +m(8); %zp
 i1vect(:, 2)=mmytry(1)*i0vect(:, 2)+mmytry(2)*i0vect(:, 1)+mmytry(3);%xp
 i1vect(:, 1)=mmytry(4) *i0vect(:, 2)+mmytry(5)*i0vect(:, 1)+mmytry(6); %yp
 i1vect(:, 3)=mmytry(7)*i0vect(:, 3) +mmytry(8); %zp
%
listpari(par+1)=pari
mytry
par
pari
mmytry
    i1vect(:, 4)= interp3(i1, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
%%%NOTE%VI = interp3(V,XI,YI,ZI) assumes X=1:N, Y=1:M, Z=1:P where [M,N,P]=size(V).

%if(par==0)
%myimage=zeros(r*c);
%myimage=i1vect(find(i1vect(:,3)==1))
%    figure; imagesc(myimage); title('interp');    
    %end

pnts=size(find(~isnan(i1vect(:, 4))), 1);
% Summing Error
err=sumskipnan((i1vect(:, 4)-i0vect(:, 4)).^2);
%errr(par-starti+1) = err/(exp(-(pnts/(2*r*c*s*1/10))^2)-1.0);
errr(par+1) = err/pnts;
end
figure; plot(listpari, errr); title(['mytry#:',int2str(mytry)]); xlabel(['m',int2str(mytry)]);
hold on
plot (mminmax(mytry, 1)    ,0,'--rs', mminmax(mytry, 2),0,'--rs')
hold on
plot (m(mytry) ,0,'-g*')

end


%errorfunction(image0,image1,m)
hold off