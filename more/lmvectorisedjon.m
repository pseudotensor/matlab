function [ m, errr, landa]=lmvectorisedjon(i0, i1, m)
% Levengberg-Marquardt Algorithm: implemented july: do not accept increment
% if the err>errp. No interpolation is used for the gradient. Instead, it
% is computed using (dI1/dx)*(dx/dm). Also the termination criterion is now
% the change in error. 
% M=lema(i0,i1,M);
%Input:
%      i0: Reference image(Gray Scale)
%      i1: Related image (Gray Scale)
%      m:  inital guess of m (column form). E.g.,  H0: m=[1 0 0 0 1 0 1 0]'
%Output:
%      M: final result: transform that takes i0 to i1.
%returns row vector of m

% jon's preperation for inverse function minimization
%i0(find(abs(i0)<.0001))=sign(i0(find(abs(i0)<.0001)))*.0001;
%i0=1/i0
%i1(find(abs(i1)<.0001))=sign(i1(find(abs(i1)<.0001)))*.0001;
%i1=1/i1
% end jon's prep
% i0=image0;
% i1=image1;
% m = [1 0 0 0 1 0 1 0]'
iteration = 0;
[r c s]=size(i0); %r=ny; c=nx; s=nz

[rp cp sp]=size(i1);
%%%! add check for the sizes to be equal


%Calculate gradient
%i1 = smooth3(i1, 'gaussian')
%i0 = smooth3(i0, 'gaussian')
%display('smoothened')

[i1gx, i1gy,  i1gz]=gradient(i1);
%[FX,FY] = gradient(F) . Fx = differences in column (from col to col) direction, fy - in row
%direction
display('gradient computed')


i1vect=zeros(rp*cp*sp, 4);

%vectorize i0
i0vect=zeros(r*c*s, 4); % vectorized i0 has 4 columns - with y, x, z coordinates and intensity value in the fourth column)
h=1:r*c*s;
i0vect(:, 1)=(mod(h-1, r)+1)';%y
i0vect(:, 2)=(floor(mod((h-1)/r, c))+1)'; %x
i0vect(:, 3)=(floor((h-1)/(c*r))+1)'; %z
i0vect(:, 4) = i0(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%what's below? revise?
%initial Landa
landa=0.0001;
pnts=size(find(~isnan((i0(:)+i1(:)))), 1);
errp=sumskipnan((i0(:)-i1(:)).^2)/pnts
errr=0
pnts=0
dm = zeros(8,1);
dontcompute=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Levengerg Marquardt algorithm
    %begin 3 cycles 
    
    m=m+dm;
    i1vect(:, 1)=m(4)*i0vect(:, 2)+m(5)*i0vect(:, 1)+m(6); %yp
    i1vect(:, 2)=m(1)*i0vect(:, 2)+m(2)*i0vect(:, 1)+m(3);%xp
    i1vect(:, 3)=m(7)*i0vect(:, 3) + m(8); %zp
    i1vect(:, 4)= interp3(i1, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');


%fid=fopen('jondata.txt','w');
%for i=1:r*c*s
%    fprintf(fid,'%g %g %g %g %g %g\n',i0vect(i,2),i0vect(i,1),i0vect(i,3),i1vect(i,2),i1vect(i,1),i1vect(i,3));
%end
%fclose(fid);
%return;


    
    A=zeros(8,8);
    b=zeros(8,1);
  
while(landa >= 1e-14 & iteration<30)

  
    %creates vectorized i1-transformed: yp, xp, zp coordinates obtained from
    %xyz. The 4th column is intensity values from i1 (interpolated)

           if(dontcompute==0)
        
                % Interpolate Gradient in Y direction
                ipgy=interp3(i1gy, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
              
				% Interpolate Gradient in X direction
				%ipgx=interp3(i1gx, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
   				ipgx=interp3(i1gx, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
                
                % Interpolate Gradient in Z direction
				%ipgz=interp3(i1gz, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
   				ipgz=interp3(i1gz, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
							
                % this is an absolutely free affine transformation, axes can not only skew at different angles,
                % but also end up parallel!
				% Compute Partial Derivate of Error
				%pem=zeros(r*c*s, 8);
        	    %pem(:, 1)=i0vect(:, 2).*ipgx; %pem(1)=j*ipgx;
				%pem(:, 2)=i0vect(:, 1).*ipgx; %pem(2)=i*ipgx;
				%pem(:, 3)=1.*ipgx; %pem(3)=1*ipgx;			
				%pem(:, 4)=i0vect(:, 2).*ipgy; %pem(4)=j*ipgy;
				%pem(:, 5)=i0vect(:, 1).*ipgy; %pem(5)=i*ipgy;
                %pem(:, 6)=1.*ipgy; %pem(6)=1*ipgy;
                %pem(:, 7)=i0vect(:, 3).*ipgz; %pem(7)=k*ipgz;
                %pem(:, 8)=1.*ipgz; %pem(8)=1*ipgz;
               
                % this forces the axes to be perpendicular.  Note that it assumes small angle changes
                %, otherwise near theta=Pi/2 we'd need to use the formula set with setting d/dm5->0 instead of d/dm4->0 as here.
                % Compute Partial Derivate of Error
				pem=zeros(r*c*s, 8);
         	    pem(:, 1)=i0vect(:, 2).*ipgx+ipgy.*(m(5)*m(2)/m(1)^2).*i0vect(:,2);
				pem(:, 2)=i0vect(:, 1).*ipgx+ipgy.*(-m(5)/m(1)).*i0vect(:,2);
				pem(:, 3)=1.*ipgx;
				pem(:, 4)=0.*ipgx; % just an easy way to get zeros
				pem(:, 5)=(i0vect(:, 1)-m(2)/m(1).*i0vect(:,2)).*ipgy;
                pem(:, 6)=1.*ipgy;
                pem(:, 7)=i0vect(:, 3).*ipgz;
                pem(:, 8)=1.*ipgz;

                % Updateing A and b
                b(1) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 1));
                b(2) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 2));
                b(3) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 3));
                b(4) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 4));
                b(5) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 5));
                b(6) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 6));
                b(7) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 7));
                b(8) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 8));
                     
                pem(any(isnan(pem)'),:) = []; %removes from PEM any rows that contain NaNs. 
                A = (pem(:, :))'*pem(:, :);
            end
               % A
               % b
            
        % Solve M's incerment
        	    dm=inv(A+landa*eye(8))*b
    mtry=m+dm;
    fprintf(1,'dm(4): %g',dm(4));
    mtry(4) = -mtry(5)*mtry(2)/mtry(1);
    fprintf(1,'mtry(4): %g',mtry(4));
    i1vect(:, 1)=mtry(4)*i0vect(:, 2)+mtry(5)*i0vect(:, 1)+mtry(6); %yp
    i1vect(:, 2)=mtry(1)*i0vect(:, 2)+mtry(2)*i0vect(:, 1)+mtry(3);%xp
    i1vect(:, 3)=mtry(7)*i0vect(:, 3) + mtry(8); %zp

%computes error between the original i0 image and i1-intensities at
%corresponding points: x&xp, y&yp, z&zp. For newm = x+ = newly updated
%pameter vector
i1vect(:, 4)= interp3(i1, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
%%%NOTE%VI = interp3(V,XI,YI,ZI) assumes X=1:N, Y=1:M, Z=1:P where [M,N,P]=size(V).
pnts=size(find(~isnan(i1vect(:, 4))), 1)
% Summing Error
err=sumskipnan((i1vect(:, 4)-i0vect(:, 4)).^2);
errr = err/pnts
    
   
    if errr>errp
          [errr errp]        
       landa=landa*10
          display('landa increased');

          dontcompute=1;
          
      else 
      [errr errp]
          landa=landa*0.1
      display('landa decreased');
   m=mtry

   errp=errr;
          dontcompute=0;   
    end
    
%
% debug!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
I3=zeros(r, c, s);
I3(:)=i1vect(:,4);
figure; imagesc(I3(:,:,1)); title(['iteration:',int2str(iteration),'dontcomp',int2str(dontcompute)]); xlabel(['error:',num2str(errr)]);

%I3(:)=ipgy(:);
%figure; imagesc(I3(:,:,1)); title(['ipgy']);
%I3(:)=ipgx(:);

%figure; imagesc(I3(:,:,1)); title(['ipgx']);
%I3(:)=ipgz(:);
%figure; imagesc(I3(:,:,1)); title(['ipgz']);

%
%
%return;
       iteration=iteration+1
end


