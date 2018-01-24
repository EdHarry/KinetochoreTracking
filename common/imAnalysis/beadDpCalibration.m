function[dp]=beadDpCalibration(image)
%beadDpCalibration calibrates the penetration depth dp from the image of a
%fluorescently coated silica bead on top of the glass; the bead is
%index-matched to the solution to prevent refraction, thus the measurement
%is performed in sucrose solution
%INPUT: image of bead
%Dinah Loerke October 11, 2004

%the following values are implemented in the function; change the
%corresponding variables to make changes (e.g. for different magnification,
%which changes the effective pixelsize, or for different wavelength)
%bead diameter: 6.84 ± 0.57 um
%camera pixel width = 6.45 um
%water r.i. =1.33
%glass coverslip r.i.= 1.523
%cell r.i. = 1.38
%wavelength = 568 nm
%assuming a magnification of 1.0x
%=>effective pixelsize of 64.5nm

pixsiz = 64.5;
%effective pixelsize
maxrad = 6840/2;
%bead radius in nm
n1 = 1.523;
% = refractive index of first medium 
n2 = 1.33;
% = refractive index of aqueous medium 
n3 = 1.36;
% = refractive index of sucrose solution (should match beads')
n4 = 1.38;
% = refractive index of cell
lambda=568;
%laser wavelanghth (in nm)
disp('current parameter values:');
disp('pixelsiz(nm), bead rad(nm), exc.wavelength (nm)');
dt=[pixsiz, maxrad, lambda];
disp(dt);
disp('n1 (glass), n2 (water), n3 (sucr), n4 (cell)');
dt=[n1, n2, n3, n4];
disp(dt);

%determine center of gravity for bead image with Gaussian fit
%first make rough estimates for starting point values of fit parameters
[imsx,imsy]=size(image);
x0i=round(imsx/2);
y0i=round(imsy/2);
offseti=min(min(image));
ai=max(max(image))-offseti;
w0i=round(x0i/4);

xw=1:imsx;
yw=1:imsy;
guess=double([ai w0i x0i y0i offseti]);
data=double(image);
distmat=double(image);
estimates = fitcurveBGauss2d(xw,yw,data,guess);

%results of Gauss fit are taken as cenetr of gravity, that is as the
%"equator" point of the bead, in reference to which all the distances are
%calculated
cg=[estimates(3) estimates(4)];
fdata=data;
%fit function and distance function are calculated
for i=1:imsx
    for j=1:imsy
        fdata(i,j)=estimates(5)+estimates(1).*exp( -(1/estimates(2))^2 * ( (i-estimates(3)).^2+(j-estimates(4)).^2 ) );
        distmat(i,j)=sqrt((i-cg(1))^2+(j-cg(2))^2);
        %if pixel lies outside of circular bead projection, distmat is set
        %to zero, marked for later point deletion
        if((distmat(i,j)*pixsiz)>maxrad)
            distmat(i,j)=0;
        end
    end
end

%original and fit function are displayed as contour plots
figure
contour(data);
hold
contour(fdata);

%axial (z-) distance is calculated from lateral distance
zdistmat=maxrad-sqrt(maxrad^2-(pixsiz.*distmat).^2);
%one-dimensional intensity and zdistance vectors are calculated from the
%two-dimensional matrices (assuming radial symmetry of the bead)
int=ones(imsx*imsy,1);
zdist=ones(imsx*imsy,1);
for i=1:imsx
    for j=1:imsy
        int((i-1)*imsy+j)=double(image(i,j));
        zdist((i-1)*imsy+j)=zdistmat(i,j);
    end
end

%remove points outside bead projection, whose zdist is set to zero
for i=(imsx*imsy):-1:1
    if (zdist(i)<1)
        zdist(i) = [];
        int(i)=[];
    end
end

%intensity versus z-distance is plotted...
figure
plot(zdist,int,'b.');
%...and fitted with an exponential function...
guess2=[(max(int)-min(int)) 5 min(int)];
estimates2 = fitcurveBexp(zdist, int, guess2);
%...which is added to the plot...
fint=estimates2(3)+estimates2(1).* exp(-zdist/estimates2(2)); 
hold
plot(zdist,fint,'r.');
xlabel('distance from coverslip (in nm)');
ylabel('fluorescence intensity (cts)');

%results: dp in sucrose
dpsuc=estimates2(2);
disp('dp in sucrose solution:');
disp(dpsuc);

%effective incidence angle corresponding to dp in sucrose solution
theta=(180/pi)*asin(sqrt(((lambda/(4*pi*dpsuc))^2+n3^2)/(n1^2)));
%dp in aequeous medium corresponding to this angle
dpaeq=lambda/(4*pi*sqrt(n1^2*(sin(pi*theta/180))^2-n2^2));
disp('dp in aequeous medium:');
disp(dpaeq);
%dp in cell (of n4=1.38) corresponding to this angle
dpcell=lambda/(4*pi*sqrt(n1^2*(sin(pi*theta/180))^2-n4^2));
disp('dp in cell:');
disp(dpcell);

title(['fitted exponential corresponds to theta=',num2str(theta),', dp=',num2str(dpaeq),' in water']);
dp=dpaeq;
end


function estimates = fitcurveBGauss2d(x,y,data,guess)
% Call fminsearch with a random starting point.
[xs,ys]=size(data);
start_point = guess;
estimates = fminsearch(@Gauss2Dfun, start_point);
% expfun accepts curve parameters as inputs and outputs sse,
% the sum of squares error for A * exp(-lambda * t) - Data.
    function sse = Gauss2Dfun(params)
        A = params(1);
        w0 = params(2);
        x0=params(3);
        y0=params(4);
        offset=params(5);
        for i=1:xs
            for j=1:ys
                FittedCurve(i,j) = offset+A .* exp( -(1/w0)^2 * ( (x(i)-x0).^2+(y(j)-y0).^2 ) );
                ErrorVector(i,j) = FittedCurve(i,j) - data(i,j);
            end
        end
        sse = sum(sum(ErrorVector .^ 2));
    end
end    


function estimates = fitcurveBexp(z, data, guess)
% Call fminsearch with a random starting point.
start_point = guess;
estimates = fminsearch(@expfun, start_point);
% expfun accepts curve parameters as inputs and outputs sse,
% the sum of squares error for A * exp(-lambda * t) - Data.
    function sse = expfun(params)
        A = params(1);
        dp = params(2);
        offset=params(3);
        FittedCurve = offset+A .* exp(-z/dp);
        ErrorVector = FittedCurve - data;
        sse = sum(ErrorVector .^ 2);
    end
end