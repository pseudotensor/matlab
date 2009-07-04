function [m, err, landa]= lmvectorised(i0, i1, m)
% Levengberg-Marquardt Algorithm
% M=lema(i0,i1,M);
%Input:
%      i0: Reference image(Gray Scale)
%      i1: Related image (Gray Scale)
%      m:  inital guess of m (column form). E.g.,  H0: m=[1 0 0 0 1 0 1 0]'
%Output:
%      M: final result: transform that takes i0 to i1.
%returns row vector of m


iteration = 1;
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


i1vect=zeros(rp*cp*sp, 4) ;

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
exit=0;
errp=0;
err=0;

while(exit==0)

	A=zeros(8,8);
	b=zeros(8,1);
	pnts=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Levengerg Marquardt algorithm
%begin 3 cycles 
%creates vectorized i1-transformed: yp, xp, zp coordinates obtained from
%xyz. The 4th column is intensity values from i1 (interpolated)
i1vect(:, 1)=m(4)*i0vect(:, 2)+m(5)*i0vect(:, 1)+m(6); %yp
i1vect(:, 2)=m(1)*i0vect(:, 2)+m(2)*i0vect(:, 1)+m(3);%xp
i1vect(:, 3)=m(7)*i0vect(:, 3) + m(8); %zp


% check this point overlapped with ip or not
			% replace by cubic later

%i1vect(:, 4)= interp3(i1, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
i1vect(:, 4)= interp3(i1, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
%%%NOTE%VI = interp3(V,XI,YI,ZI) assumes X=1:N, Y=1:M, Z=1:P where [M,N,P]=size(V).

pnts=size(find(~isnan(i1vect(:, 4))), 1);

           % Summing Error
err=sumskipnan((i1vect(:, 4)-i0vect(:, 4)).^2);
              
                
				% Interpolate Gradient in Y direction
				%ipgy=interp3(i1gy, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
                ipgy=interp3(i1gy, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
              
				% Interpolate Gradient in X direction
				%ipgx=interp3(i1gx, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
   				ipgx=interp3(i1gx, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
                
           
                % Interpolate Gradient in Z direction
				%ipgz=interp3(i1gz, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
   				ipgz=interp3(i1gz, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
				
								
				% Compute Partial Derivate of Error
				pem=zeros(r*c*s, 8);
         	    pem(:, 1)=i0vect(:, 2).*ipgx; %pem(1)=j*ipgx;
				pem(:, 2)=i0vect(:, 1).*ipgx; %pem(2)=i*ipgx;
				pem(:, 3)=1.*ipgx; %pem(3)=1*ipgx;			
				pem(:, 4)=i0vect(:, 2).*ipgy; %pem(4)=j*ipgy;
				pem(:, 5)=i0vect(:, 1).*ipgy; %pem(5)=i*ipgy;
                pem(:, 6)=1.*ipgy; %pem(6)=1*ipgy;
                pem(:, 7)=i0vect(:, 3).*ipgz; %pem(7)=k*ipgz;
                pem(:, 8)=1.*ipgz; %pem(8)=1*ipgz;
                pem;

% Updateing A and b
                b(1) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 1));
                b(2) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 2));
                b(3) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 3));
                b(4) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 4));
                b(5) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 5));
                b(6) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 6));
                b(7) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 7));
                b(8) = -sumskipnan((i1vect(:, 4)-i0vect(:, 4)).*pem(:, 8));
                 b;
                 
                pem(any(isnan(pem)'),:) = []; %removes from PEM any rows that contain NaNs. 
                A = (pem(:, :))'*pem(:, :);
                A;


           
      
%%% end 3 cycles

A;
pnts;
err;
err = err/pnts;
%err=err/(r*c*s)
	% Change landa

	if err>errp
		landa=10*landa
	else
		landa=0.1*landa
	end         

	% Check Stop criterion
    
	if (landa <= 1e-14) | (iteration>100)
	       m = m' % returns row- vector of m
 
           exit=1;      
	else
	    % Solve M's incerment
	    dm=inv(A+landa*eye(8))*b;
	    % Update M
	    m=m+dm %both m and dm are columns
	    landa 
	    [err errp]
        iteration=iteration+1;
	    iteration
	    errp=err;
	    err=0;
    end
end
fid=fopen('C:\temp\test.txt', 'w')
fprintf(fid, '%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f' , m, err, landa)
fclose(fid)