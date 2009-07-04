%function [I2,I3]=myNewaffine(I0, m)
function [I2]=myNewaffine(I0, m)
%My own affine transformation routine that preserves the size of the image
%minn is the background value (put ANYTHING). Now performs 3D affine!

%pnts=0;
%[ilength jlength klength] = size(I0);
%I2 = zeros(ilength, jlength, klength);
%for k = 1:klength
%for i = 1:ilength
%	for j = 1:jlength
	
		%compute the transformed coordinates using the current estimates of the parameters.
 %         	x=-(m(3)*m(5)-m(2)*m(6)-m(5)*j+m(2)*i)/(-m(2)*m(4)+m(1)*m(5))
  %        	y=-(m(3)*m(4)-m(1)*m(6)-m(4)*j+m(1)*i)/(m(2)*m(4)-m(1)*m(5))
   %         z = k; %/m(7)
            	% check whether this new point is still within image frameworks
		%if ((x>=1) & (x<= jlength) & (y>=1) & (y<= ilength)&(z>=1) & (z<= klength))
   %        pnts=pnts+1;
          
             %Interpolate Intensity
		%	dx=x-floor(x);
		%		dy=y-floor(y);
		%		a0=I0(floor(y),floor(x));
		%		a1=I0(floor(y),ceil(x));
		%		a2=I0(ceil(y),floor(x));
			%	a3=I0(ceil(y),ceil(x));
		%		a01=a0+(a1-a0)*dx;
	%			a23=a2+(a3-a2)*dx;
%				intensity=a01+(a23-a01)*dy;
   %        I2(i, j, k) = interp3(I0, x, y, z, 'linear');      
   %end
   %end
   %end
   %end
[r c s] = size(I0);
i0vect=zeros(r*c*s, 3);
I2=zeros(r, c, s);

i2vect=zeros(r*c*s, 4); % vectorized i0 has 4 columns - with y, x, z coordinates and intensity value in the fourth column)
h=1:r*c*s;
i2vect(:, 2)=(floor(mod((h-1)/r, c))+1)'; %xp
i2vect(:, 1)=(mod(h-1, r)+1)';%yp
i2vect(:, 3)=(floor((h-1)/(c*r))+1)'; %zp

% coordinates transforming x',y',z' to x,y,z where we assume I0 is located initially in x',y',z'.
i0vect(:, 2)=-(m(3)*m(5)-m(2)*m(6)-m(5)*i2vect(:, 2)+m(2)*i2vect(:, 1))./(-m(2)*m(4)+m(1)*m(5));%x
i0vect(:, 1)=-(m(3)*m(4)-m(1)*m(6)-m(4)*i2vect(:, 2)+m(1)*i2vect(:, 1))./(m(2)*m(4)-m(1)*m(5));%y
i0vect(:, 3) = (i2vect(:, 3) -m(8))./m(7);%z



%fid=fopen('jondata.txt','w');
%for i=1:r*c*s
%    fprintf(fid,'%g %g %g %g %g %g\n',i0vect(i,2),i0vect(i,1),i0vect(i,3),i2vect(i,2),i2vect(i,1),i2vect(i,3));
%end
%fclose(fid);

% this gives new image in x,y,z coming from x',y',z'
I2(:)= interp3(I0, i0vect(:, 2), i0vect(:, 1), i0vect(:, 3), 'linear');

% jon stuff
% now let's test our back transform.  Assume now i2vect is x,y,z
% and apply transform from I2(x,y,z) -> I3(x',y',z')
i1vect(:, 1)=m(4)*i2vect(:, 2)+m(5)*i2vect(:, 1)+m(6); %yp
i1vect(:, 2)=m(1)*i2vect(:, 2)+m(2)*i2vect(:, 1)+m(3);%xp
i1vect(:, 3)=m(7)*i2vect(:, 3) + m(8); %zp
%
I3=zeros(r, c, s);
I3(:)= interp3(I2, i1vect(:, 2), i1vect(:, 1), i1vect(:, 3), 'linear');
% now I0(x,y,z) and I3(x',y',z') should be comparable in viewing images on each one's coordinates
figure; imagesc(I0(:,:,1)); title('0');
figure; imagesc(I2(:,:,1)); title('2');
figure; imagesc(I3(:,:,1)); title('3');
