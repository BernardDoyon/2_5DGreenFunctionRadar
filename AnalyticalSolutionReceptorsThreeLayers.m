% ======================================================================
% Copyright (c) June 2019, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================

function AnalyticalSolutionReceptorsThreeLayers(mu_r,epsilon_r, ...
                                        sigma,delta_z,z0,sx,sz,f,q,...
                                        rx,rz,ry,pola_source,FileName)

%   ANALYTICAL SOLUTION FOR THE ELECTRIC FIELD FOR
%       THE THREE LAYER MODEL
%
%      We use the convention for the time oscillating part
%      E = E_0*exp(-i*omega*t)
%
%      Solution is coded only for the midle layer
%   
%   INPUT:
%
%   mu_r: Relative permeability vector (from top to bottom, only different
%                                       layers are possible)
%   epsilon_r: Relative permittivity vector
%   sigma: Electric conductivity vector (Siemens/meter)
%   delta_z: grid spacing in z
%   z0: Top left corner z coordinate of dommain
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

% Top layer  (index 2)
% Middle layer (index 1)
% Bottom layer (index 3)

NbSources = length(sx);
NbReceptors = length(rx);

indexInterface = find(diff(epsilon_r));

epsilon2 = epsilon_r(1)*epsilon0;
% z1 = z coordonnate of beginning of middle layer
z1 = indexInterface(1)*delta_z + z0;

b = (indexInterface(2)-indexInterface(1))*delta_z;
epsilon1 = epsilon_r(indexInterface(1)+1)*epsilon0;
epsilon3 = epsilon_r(end)*epsilon0;
z2= b+z1; % z2: z coordonnate of beginning of bottom layer

flagReceptor = zeros(1,NbReceptors); % flag to inidcate in which layer
% the receptors are
for i = 1:NbReceptors
    if (rz(i)<z1)
        flagReceptor(i) = 2;
    elseif (rz(i)<z2)
        flagReceptor(i) = 1;
    else
        flagReceptor(i) = 3;
    end
end

mu2 = mu_r(1)*mu0;
mu1 = mu_r(round(length(mu_r)/2))*mu0;
mu3 = mu_r(end)*mu0;



sigma2 = sigma(1);
sigma1 = sigma(round(length(epsilon_r)/2));
sigma3 = sigma(end);


Nf = length(f(:,1));
omega = 2*pi*(f(:,1) + 1i*q)*1e6;


EzR = zeros(Nf,length(rx));
ExR = zeros(Nf,length(rx));
EyR = zeros(Nf,length(rx));

if (pola_source == 1) % the solution is coded only for source polarized 
                           % z direction
     for w =1:Nf
         
         % For every complexe frequency, we calculate the E field
         % We use the convention for the time oscillating part
         %  E = E_0*exp(-i*omega*t)
         %
         % This impose the wave number the wave to propagate in space with
         % E = E_0 exp(i*k*r) (the wave number k is POSITIVE)
         
         YAdm2 = -1i*omega(w)*epsilon2 + sigma2;
         ZImp2 = -1i*omega(w)*mu2;
         k2 = abs(real(sqrt(-YAdm2*ZImp2))) + 1i*abs(imag(sqrt(-YAdm2*ZImp2)));
       
         YAdm1 = -1i*omega(w)*epsilon1 + sigma1;
         ZImp1 = -1i*omega(w)*mu1;
         k1 = abs(real(sqrt(-YAdm1*ZImp1))) + 1i*abs(imag(sqrt(-YAdm1*ZImp1)));
         
         YAdm3 = -1i*omega(w)*epsilon3 + sigma3;
         ZImp3 = -1i*omega(w)*mu3;
         k3 = abs(real(sqrt(-YAdm3*ZImp3))) + 1i*abs(imag(sqrt(-YAdm3*ZImp3)));
         
         
         YAdm = [YAdm1 YAdm2 YAdm3];
         ZImp = [ZImp1 ZImp2 YAdm3];
       
         Kc = max(real(sqrt(-YAdm.*ZImp)));
         
         %KrMax = pi/delta_z;
         KrMax = 10*Kc;
         
         Nk = 1000;
         Kr = KrMax*linspace(0,1,Nk+1); % kr de 0 à kr_max
         
         M = 1/4/pi/YAdm1;
         
         for r = 1:NbReceptors
             
             switch flagReceptor(r)
                 case 1
                     for n = 1:NbSources
                         Sum1 = 0;
                         h = sz(n) - z1; %z distance between source and
                         % first interface (layer2-layer1)
                         
                         X = (rx(r) - sx(n));
                         Z = (rz(r) - sz(n)) + h;
                         Y = ry;    % The sources are located at y = 0
                         rho = sqrt(X.*X + Y*Y);
                         R = sqrt(rho.*rho + (Z-h).*(Z-h));
                         delta_k = Kr(2)-Kr(1);
                         
                         for k = 2:Nk
                             k1z = sqrt(k1*k1 - Kr(k)*Kr(k)); % tau1 Tyras
                             k2z = sqrt(k2*k2 - Kr(k)*Kr(k)); % tau2 Tyras
                             k3z = sqrt(k3*k3 - Kr(k)*Kr(k)); % tau3 Tyras
                             R12 = ((k2/k1)^2*k1z - k2z)/((k2/k1)^2*k1z + k2z);
                             R13 = ((k3/k1)^2*k1z - k3z)/((k3/k1)^2*k1z + k3z);
                             A = R12*(1+R13*exp(2i*k1z*(b-h)))/(1-R12*R13*exp(2i*k1z*b))*exp(1i*k1z*h)/k1z;
                             B = R13*exp(2i*k1z*b)*(1+R12*exp(2i*k1z*h))/(1-R13*R12*exp(2i*k1z*b))*exp(-1i*k1z*h)/k1z;
                             Fc = (A*exp(1i*k1z*Z) + B*exp(-1i*k1z*Z));
                             d_sum = Fc.*besselj(0,Kr(k)*rho)*(Kr(k))^3*delta_k;
                             
                             Max_dsum = abs(d_sum);
                             Sum1 = Sum1 + d_sum;
                             
                             if (((Max_dsum/abs(Sum1) < 1e-12)||(Max_dsum < 1e-12))&&(Kr(k) > Kc))
                                 break;
                             end
                         end
                         EzR(w,r) =  EzR(w,r) + 1i*M*Sum1 + M*exp(1i*k1*R).*(3*(Z-h).*(Z-h)./R.^5 - 3i*k1*(Z-h).*(Z-h)./R.^4 ...
                                         - (k1^2*(Z-h).*(Z-h) + 1)./R.^3 + 1i*k1./R.^2 + k1^2./R);
                     end
                 case 2
                     %{
                      for n = 1:NbSources
                          Sum1 = 0;
                          h = sz(n) - z1; %z distance between source and
                          % first interface (layer2-layer1)
                          
                          X = (rx(r) - sx(n));
                          Z = (rz(r) - sz(n)) + h;
                          Y = ry;    % The sources are located at y = 0
                          rho = sqrt(X.*X + Y*Y);
                          delta_k = Kr(2)-Kr(1);
                          
                          for k = 2:Nk
                              k1z = sqrt(k1*k1 - Kr(k)*Kr(k)); % tau1 Tyras
                              k2z = sqrt(k2*k2 - Kr(k)*Kr(k)); % tau2 Tyras
                              k3z = sqrt(k3*k3 - Kr(k)*Kr(k)); % tau3 Tyras
                              R12 = ((k2/k1)^2*k1z - k2z)/((k2/k1)^2*k1z + k2z);
                              R13 = ((k3/k1)^2*k1z - k3z)/((k3/k1)^2*k1z + k3z);
                              T12 = 2*exp(i*k1z*h)/((k2/k1)^2*k1z + k2z)*((1+R13*exp(2*i*k1z*(b-h)))/(1-R12*R13*exp(2*i*k1z*b)));
                          
                              d_sum = T12*exp(-1i*k2z*Z)*besselj(0,Kr(k)*rho)*(Kr(k))^3*delta_k;
                              
                              Max_dsum = abs(d_sum);
                              Sum1 = Sum1 + d_sum;
                              
                              if (((Max_dsum/abs(Sum1) < 1e-12)||(Max_dsum < 1e-12))&&(Kr(k) > Kc))
                                  break;
                              end
                          end
                          EzR(w,r) =  EzR(w,r) + 1i*M*Sum1;
                      end
                     %}
                      EzR(w,r)=NaN;   % solution not available in this layer
                 case 3
                     %{
                     for n = 1:NbSources
                         Sum1 = 0;
                         h = sz(n) - z1; %z distance between source and
                         % first interface (layer2-layer1)
             
                         X = (rx(r) - sx(n));
                         Z = (rz(r) - sz(n)) + h;
                         Y = ry;    % The sources are located at y = 0
                         rho = sqrt(X.*X + Y*Y);
                         delta_k = Kr(2)-Kr(1);
                         
                         for k = 2:Nk
                             k1z = sqrt(k1*k1 - Kr(k)*Kr(k)); % tau1 Tyras
                             k2z = sqrt(k2*k2 - Kr(k)*Kr(k)); % tau2 Tyras
                             k3z = sqrt(k3*k3 - Kr(k)*Kr(k)); % tau3 Tyras
                             R12 = ((k2/k1)^2*k1z - k2z)/((k2/k1)^2*k1z + k2z);
                             R13 = ((k3/k1)^2*k1z - k3z)/((k3/k1)^2*k1z + k3z);
                             T13 = 2*exp(i*k1z*(b-h))/((k3/k1)^2*k1z + k3z)*((1+R12*exp(2*i*k1z*h))/(1-R12*R13*exp(2*i*k1z*b)));
                             d_sum = T13*exp(1i*k3z*(Z-b))*besselj(0,Kr(k)*rho)*(Kr(k))^3*delta_k;
                             
                             Max_dsum = abs(d_sum);
                             Sum1 = Sum1 + d_sum;
                             
                             if (((Max_dsum/abs(Sum1) < 1e-12)||(Max_dsum < 1e-12))&&(Kr(k) > Kc))
                                 break;
                             end
                         end
                         EzR(w,r) =  EzR(w,r) + 1i*M*Sum1;
                     end
                     %}
                     EzR(w,r)=NaN;   % solution not available in this layer
             end
         end
         ExR(w,:) = NaN; % Ex is not coded. 
         EyR(w,:) = NaN;
     end
else
    EzR(:,:) = NaN;  %#ok<*NASGU>
    ExR(:,:) = NaN;
    EyR(:,:) = NaN;
end
save(FileName,'EzR','ExR','EyR','f','q','sx','sz','rx','rz','ry');
end
 
  



            
          

           

