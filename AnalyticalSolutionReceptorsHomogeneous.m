% ======================================================================
% Copyright (c) June 2019, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================

function AnalyticalSolutionReceptorsHomogeneous(mu_r,epsilon_r, sigma,...
                        sx,sz,f,q,rx,rz,ry,pola_source,FileName)

%   ANALYTICAL SOLUTION FOR THE ELECTRIC FIELD WITH
%       CONSTANT MU, EPSILON AND SIGMA 
%
%      We use the convention for the time oscillating part
%      E = E_0*exp(-i*omega*t) 
%
%   INPUT:
%
%   mu_r: Relative permeability 
%   epsilon_r: Relative permittivity 
%   sigma: Electric conductivity (Siemens/meter)
%   sx: x position of sources (meters)
%   sz: z position of sources (meters)
%       Note: 
%           If N sources than sx(1:N) and sz(1:N)
%           (sources are loacated at y = 0)
%   frequencies: frequency vector of length Nf (MHz)
%   q: Imaginary component of frequency (MHz)  
%   rx: x positions of receptors (meters)
%   rz: z positions of receptors (meters)
%   ry: y positions of receptors (meters)
%       Note:
%           If N receptors than srx(1:N) and rz(1:N)
%   pola_source: (Choices are: 0, 1 or 2)
%       Note:
%           - pola_source = 0 => Source Polarized in the "y" direction
%           - pola_source = 1 => Source Polarized in the "z" direction
%           - pola_source = 2 => Source Polarized in the "x" direction
%   Filename: Filename to store the solution (.mat file) 
%
%   OUTPUT (in the Filename.mat file)
%
%   ExR(Nf,Number of receptors): x component of the E field
%   EzR(Nf,Number of receptors): z component of the E field
%   EyR(Nf,Number of receptors): y component of the E field
%
%   Example:
%           ExR(3,2) is the complex x componants of the Green function at
%           receptor index 2
%

     mu0 = 4*pi*1e-7;
     epsilon0 = 1/36/pi*1e-9;
     
     epsilon = epsilon_r*epsilon0;
     mu = mu_r*mu0;

     Nf = length(f(:,1));

     omega = 2*pi*(f(:,1) + 1i*q)*1e6;
     
     ExR = zeros(Nf,length(rx));
     EzR = zeros(Nf,length(rx));
     EyR = zeros(Nf,length(rx));
     
     NbSources = length(sx);
     
     for w =1:Nf
       
         for n = 1:NbSources
             X = (rx - sx(n));
             Z = (rz - sz(n));
             Y = ry;    % The sources are located at y = 0
             
             % For every complexe frequency, we calculate the E field
             % We use the convention for the time oscillating part
             %  E = E_0*exp(-i*omega*t)
             %
             % This impose the space propagation to be:
             % E = E_0 exp(i*k*r) (the wave number k is POSITIVE)
             
             YAdm = -1i*omega(w)*epsilon + sigma;
             ZImp = -1i*omega(w)*mu;
             k = abs(real(sqrt(-YAdm*ZImp))) + 1i*abs(imag(sqrt(-YAdm*ZImp)));
             
             R = sqrt(X.*X + Z.*Z + Y*Y); % Distance from source
             
             % Convention E = E_0*exp(-i*omega*t)
             % => E = E_0 exp(i*k*r) (the wave number k is POSITIVE)
             
             A =  -R.^2*k^2 + 3 - 3i*k*R;
             B =  R.^2*k^2 - 1 + 1i*k*R;
             
             switch(pola_source)
                 case 1
                     ExR(w,:) = ExR(w,:) + (X.*Z./(R.*R).*A).*exp(1i*k*R)/4/pi/YAdm./R.^3;
                     EyR(w,:) = EyR(w,:) + (Y.*Z./(R.*R).*A).*exp(1i*k*R)/4/pi/YAdm./R.^3;
                     EzR(w,:) = EzR(w,:) + (Z.*Z./(R.*R).*A + B).*exp(1i*k*R)/4/pi/YAdm./R.^3;
                 case 2
                     EzR(w,:) = EzR(w,:) + (X.*Z./(R.*R).*A).*exp(1i*k*R)/4/pi/YAdm./R.^3;
                     ExR(w,:) = ExR(w,:) + (X.*X./(R.*R).*A + B).*exp(1i*k*R)/4/pi/YAdm./R.^3;
                     EyR(w,:) = EyR(w,:) + (Y.*X./(R.*R).*A).*exp(1i*k*R)/4/pi/YAdm./R.^3;
                 case 0
                     ExR(w,:)= ExR(w,:) + (X.*Y./(R.*R).*A).*exp(1i*k*R)/4/pi/YAdm./R.^3;
                     EzR(w,:) = EzR(w,:) + (Y.*Z./(R.*R).*A).*exp(1i*k*R)/4/pi/YAdm./R.^3;
                     EyR(w,:) = EyR(w,:) + (Y.*Y./(R.*R).*A + B).*exp(1i*k*R)/4/pi/YAdm./R.^3;
             end
             % If one wants to have the convention E = E_0*exp(i*omega*t)
             % => E = E_0 exp(-i*k*r) (with wave number k POSITIVE)
             
             % A =  -R.^2*k^2 + 3 + 3i*k*R;
             % B =  R.^2*k^2 - 1 - 1i*k*R;
             % Ex = -R2omega_p(w)*(X.*Z./(R.*R).*A).*exp(-1i*k*R)/4/pi/YAdm./R.^3;
             % Ez = -R2omega_p(w)*(Z.*Z./(R.*R).*A + B).*exp(-1i*k*R)/4/pi/YAdm./R.^3;
         end
         for j = 1:length(rx) 
             if (abs(ExR(w,j))<1e-15)
                 ExR(w,j) = NaN;
             end
             if (abs(EzR(w,j))<1e-15)
                 EzR(w,j) = NaN;
             end
             if (abs(EyR(w,j))<1e-15)
                 EyR(w,j) = NaN;
             end
         end
     end
     save(FileName,'ExR', 'EzR','EyR','f','q','sx','sz','rx','rz','ry');
end
 
    