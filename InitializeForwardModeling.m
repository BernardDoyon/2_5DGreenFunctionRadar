function InitializeForwardModeling(FileNameModel,Nf,fmin,fmax,idexFrequency,q,lx,lz,...
                           sx,sz,rx,rz,ry,lpmlx,lpmlz,pola_source,...
                           pola_receptor,delta_x,delta_z,NbGridPointPerWavelength,...
                           NbGridSizes,InterpolationMethod,flagPlotModel) %#ok<INUSL>

%   This function will create a FileNameModelForwInfo.mat file in the 
%   CurrentFolder\Model directory. This .mat file will contain all the 
%   information needed to calculate the numerical solution for the Green 
%   function with the help of EM_NumericalGreenFunction.m file
%   (see main.m file)
%
%   The user has to provide some informations about the model, informations 
%   listed below:
%
%   1) FileNameModel: 
%
%       FileNameModel is a .mat file with the 3 matrices
%       defining the model. The 3 matrices are:
%
%       1) epsilonMatrix(NzCell,NxCell)
%       2) muMatrix(NzCell,NxCell)
%       3) sigmaMatrix(NzCell,NxCell)
%
%       with
%           NzCell: Number of cell in the z direction
%           NxCell: Number of cell in the x direction
%
%       All the other parameters of the model (for instance the length
%       of the domain or the position of the sources) are specified
%       with the
%
%   2) Nf:  Number of positive frequencies to describe the initial 
%           time domain source pulse 
%    
%   3) fmax: Maximum value of the frequency (MHz) 
%
%   4) fmin: Minimum value of the frequency (MHz)
%
%   5) idexFrequency 
%                   Index of the frequencies for which a Green function has
%                   to be calculated.
%                   
%                   Example:
%
%                   If Nf = 10; (10 positive frequencies) and fmax = 100 
%                   and fmin = 0, the positive frequencies describing the
%                   initial time domain trace of the source would be
%                   f =[10 20 30 40 50 60 70 80 90 100] MHz
%                   If we set idexFrequency = [3 6 9], the Green function
%                   will be calculated for f = 30, 60 and 90 MHz. 
%        
%   6) q:    Imaginary component of frequency (MHz)  
%   7) lx:   x length of domain (meters)
%   8) lz:   z length of domain (meters)
%   9) sx:   x position of sources (meters)
%   10) sz:   z position of sources (meters)
%   11) rx:   x position of receptors (meters)
%   12) rz:   z position of receptors (meters)
%   13) ry:   y position of receptors (meters)
%   14) lpmlx: x length of pml region (meters)
%   15) lpmly: y length of pml region (meters)
%   16) pola_source: (Choices are: 0, 1 or 2)
%       Note:
%           - pola_source = 0 => Source Polarized in the "y" direction
%           - pola_source = 1 => Source Polarized in the "z" direction
%           - pola_source = 2 => Source Polarized in the "x" direction
%
%   17) pola_receptor (Choices are: 0, 1, 2 or 3)
%       Note:
%           The numerical solution is obtained by summation. The summation 
%           will stop when the relative correction term is smaller than 1% 
%           for the  
%           - "y" component of the Green function (pola_receptor = 0), 
%           - "z" component of the Green function (pola_receptor = 1), 
%           - "x" component of the Green function (pola_receptor = 2), 
%           - all components of the Green function (pola_receptor = 3).
%  
%   18) Delta: Grid spacing (meters)
%       Note:
%           If Delta = 0, a grid spacing will be calculated in order that 
%           the minimum number of points per wavelength requested by the
%           user is respected (see variable "NbGridPointPerWavelength" 
%           below). For Delta = 0, grid spacing can be different from one 
%           frequency to the other. A value of Delta is first calculated 
%           for the highest frequency. The value of Delta is double every 
%           time the criteria for the NbGridPointPerWavelength is 
%           respected. It is possible to restrict the number of different 
%           Delta values (see variable NbGridSizes below
%
% 19) NbGridPointPerWavelength: The minimum number of points per wavelength
%                               to be respected for a given frequency.
%   
% 20) NbGridSizes: The maximum nuber of different Delta values. 
%
% 21) InterpolationMethod:  The interpolation method
%
% 22)* flagPlotModel: When this flag is set to 1, a color image of the
%                   initial model for epsilon and sigma are plot for the
%                   initial domain. Color images of the interpolated
%                   models for epsilon and sigma are also plot for each
%                   numerical grid.                

FileNameOut = strcat(FileNameModel,'InfoForwMod');

load(FileNameModel,'epsilonMatrix','muMatrix','sigmaMatrix');

% Evaluate the grid size for every frequency
% This grid size will be used to calculate
% the Numerical Green's Functions by finite differences 
f = ones(Nf,3);
temp = linspace(fmin,fmax,Nf+1);
f(:,1) = temp(2:Nf+1); 
f(:,3) = 0;
f(idexFrequency,3)=1;   % f(:,3) is a flag to indicate if the Green 
                        % Function has to be evaluated for this frequency

if ((delta_x==0)&&(delta_z == 0)) % if no Delta is specified
% Calculate Grid size to respect the number of points per wavelength                 
     [delta_x,delta_z, f] = GridSpacing(epsilonMatrix,muMatrix,sigmaMatrix,...
                        f,q,lx,lz,NbGridPointPerWavelength,NbGridSizes);
else
    [Nz, Nx]=size(epsilonMatrix);
    delta_xMod = lx/Nx;
    delta_zMod = lz/Nz;
    if(delta_x>delta_xMod)
        delta_x = delta_xMod;
    end
    if(delta_z>delta_zMod)
        delta_z = delta_zMod;
    end
end

%%%  Adjust the domain values to the selected grid %%%%%%%%
%   
% sx: x coordonnates of the sources (in meters) specified by the user
% sz: z coordonnates of the sources (in meters) specified by the user             
% sxGrid coordonnates of the sources adjusted to the grid size 
% szGrid coordonnates of the sources adjusted to the grid size
%   Note: source position is adjusted:
%       - on a grid node for source polarized in the z direction,
%       - in the center of a cell for source polarized in the x direction,
%       - on the edge of a cell for source polarized in the y direction
%
% lx:  x  length (in meters) of the modeling domain specified by the user
% lz:  z  length (in meters) of the modeling domain specified by the user
% lx_Grid: x adjusted length of the domain
% lz_Grid: z adjusted length of the domain
%   Note Domain size is adjusted so;
%       - lx_Grid/delta_x = integer for all Delta values 
%       - lz_Grid/delta_z = integer for all Delta values
    
[sx_idex,sz_idex,sx_Grid,sz_Grid, x0, z0,lx_Grid,lz_Grid,Npmlx,Npmlz] = ...
                                               AdjustDomain(lx,lz,sx,sz,...
                                               lpmlx,lpmlz,delta_x,delta_z,...
                                               pola_source); %#ok<*ASGLU>

% Interpolate with the initial parameter model to get
% parameter values on a grid
[epsilonMatrix_Grid, sigmaMatrix_Grid, muMatrix_Grid, Nx, Nz]  = ...
                        Interpol(lx,lz,lx_Grid,lz_Grid,x0,z0,delta_x,...
                                    delta_z,epsilonMatrix,muMatrix,...
                                    sigmaMatrix,InterpolationMethod);
                                
% We plot the parameter model and show the position of the sources for
% every grid sizes if flagPlotModel = 1;

fprintf('plot \n');
if (flagPlotModel) 
    PlotModelAndGrid(lx,lz,epsilonMatrix,sigmaMatrix,sx,sz,rx,rz,...
                        lx_Grid,lz_Grid,Nx,Nz,epsilonMatrix_Grid,...
                         sigmaMatrix_Grid,sx_Grid,sz_Grid,delta_x,...
                         x0,z0,f,Npmlx,Npmlz)
end

fprintf('save \n');
save(FileNameOut, 'f','q','pola_source',...
                   'pola_receptor','delta_x','delta_z', 'Nx', 'Nz','x0',...
                    'z0','Npmlx','Npmlz', 'sx_idex', 'sz_idex','rx',...
                    'rz','ry', 'epsilonMatrix_Grid',...
                    'sigmaMatrix_Grid','muMatrix_Grid');
end



function [delta_x, delta_z, f] = GridSpacing(epsilonMatrix,muMatrix,sigmaMatrix,...
                           f,q,lx,lz,NbGridPointPerWavelength,NbGridSizes)

  %%%%%%  
  % A grid size is evaluated for each frequency to respect the minimum 
  % number of gridpoints per wavelength. The highest value of the frequency 
  % will set the smallest value of grid size. The grid size will be
  % increased by a factor of 2 whenever the minimum number of gridpoints 
  % per wavelength is respected (for a given frequency). 
  % We limit the number of different grid size to 3. 
  % 
  %  
  %
    
  Kmin = 1/NbGridPointPerWavelength;
  mu0 = 4*pi*1e-7;
  epsilon0 = 1/36/pi*1e-9;
  mu = muMatrix*mu0;
  epsilon = epsilonMatrix*epsilon0;
  sigma = sigmaMatrix;
  
  Nf = length(f(:,1));
  
  [Nz, Nx]=size(epsilonMatrix);
  
  delta_xMod = lx/Nx;
  
  delta_zMod = lz/Nz;
  
  omega = 2*pi*f(:,1)*1e6 + 1i*q;
  Delta = zeros(1,Nf);
  delta_x = zeros(1,Nf);
  delta_z = zeros(1,Nf);
  for w = 1:Nf
      YAdm = -1i*omega(w)*epsilon + sigma;
      ZImp = -1i*omega(w)*mu;
      beta = abs(real(sqrt(-YAdm.*ZImp)));
      Delta(1,w) = round(2*pi*Kmin/max(max(beta)),1,'significant');
  end
  
  Delta_min = min(Delta);  % smallest value of delta (high frequency)
  % to respect the number of points per
  % wavelength
  
  
  % set values for delta_x
  temp = delta_xMod;
  while (temp > Delta_min)
      temp = temp/2;
  end
  Delta_min = temp;
  
  if (2*Delta_min>delta_xMod)
      DeltaNext_min = Delta_min;
  else
      DeltaNext_min = 2*Delta_min;
  end
     
  idexChangeGrid = 0;
  idexGrid_x = zeros(1,Nf);
  for w = Nf:-1:1
      if ((Delta(1,w)>= DeltaNext_min)&&(idexChangeGrid<NbGridSizes-1))
          delta_x(1,w) = DeltaNext_min;
          Delta_min = DeltaNext_min;
          if (2*Delta_min>delta_xMod)
              DeltaNext_min = Delta_min;
          else
              DeltaNext_min = 2*Delta_min;
          end
          idexChangeGrid = idexChangeGrid + 1;
         % f(w,2) = NbGridSizes - idex+1;
          idexGrid_x(w) = NbGridSizes - idexChangeGrid;
      else
          delta_x(1,w) = Delta_min;
          %f(w,2) = NbGridSizes - idex+1;
          idexGrid_x(w) = NbGridSizes - idexChangeGrid;
      end
  end
  
   % set values for delta_z
  Delta_min = min(Delta);
  temp = delta_zMod;
  while (temp > Delta_min)
      temp = temp/2;
  end
  Delta_min = temp;
  
  if (2*Delta_min>delta_zMod)
      DeltaNext_min = Delta_min;
  else
      DeltaNext_min = 2*Delta_min;
  end
 
  idexChangeGrid = 0;
  idexGrid_z = zeros(1,Nf);
  for w = Nf:-1:1
      if ((Delta(1,w)>= DeltaNext_min)&&(idexChangeGrid<NbGridSizes-1))
          delta_z(1,w) = DeltaNext_min;
          Delta_min = DeltaNext_min;
          if (2*Delta_min>delta_zMod)
              DeltaNext_min = Delta_min;
          else
              DeltaNext_min = 2*Delta_min;
          end
          idexChangeGrid = idexChangeGrid + 1;
          idexGrid_z(w) = NbGridSizes - idexChangeGrid;
      else
          delta_z(1,w) = Delta_min;
          idexGrid_z(w) = NbGridSizes - idexChangeGrid;
      end
  end
  % Define the different grid sizes for each frequency
  % For each frequency, we use an index to indicate the grid to use
  % This index is located in f(:,2);
  temp_x = delta_x;
  temp_z = delta_z;
  idexGrid = 1;
  delta_x = delta_x(1);
  delta_z = delta_z(1);
  f(1,2) = idexGrid;
 
  for w = 2:Nf
      if (~(idexGrid_x(w)==idexGrid_x(w-1))||~(idexGrid_z(w)==idexGrid_z(w-1)))
          idexGrid = idexGrid+1;
          delta_x(idexGrid) = temp_x(w);
          delta_z(idexGrid) = temp_z(w);
          f(w,2) = idexGrid;
      else
          f(w,2) = idexGrid;
      end
  end
  
end


function [sx_idex,sz_idex,sxGrid,szGrid, x0, z0,lx_Grid, lz_Grid, ...
            Npmlx,Npmlz]  = AdjustDomain(lx,lz,sx,sz,lpmlx,lpmlz,...
                                            delta_x,delta_z,pola_source)
     % Try to choose a grid size so the sources are located:
     %  - On a gride node for sources polarized in the z direction
     %  - In the center of a cell for sources polarized in the x direction
     %  - On the edge of a cell for sources polarized in the y direction
    
    
     NbGrid = length(delta_x);

     delta_x_max = max(delta_x);
     delta_z_max = max(delta_z);
     
     NbSources = length(sx);
     sx_idex = zeros(NbGrid,NbSources);
     sz_idex = zeros(NbGrid,NbSources);

     Npmlx = round(lpmlx/delta_x_max)*delta_x_max./delta_x;
     Npmlz = round(lpmlz/delta_z_max)*delta_z_max./delta_z;
     
     x0 = zeros(1,NbGrid)-Npmlx.*delta_x;
     z0 = zeros(1,NbGrid)-Npmlz.*delta_z; 
     
     % Adjust length of domain so that lx_Grid = Nx*delta_x
     % with Nx an integer. Same for lz_Grid...
     lx_Grid = round(lx./delta_x_max).*delta_x_max + 2*abs(x0(1));
     lz_Grid = round(lz./delta_z_max).*delta_z_max + 2*abs(z0(1));
     
     % Adjust position of sources:
     % (sxGrid and szGrid variables bellow)
     %
     % - In a corner of a cell if source polarized in z direction  
     % - In the middle of the edge of a cell if source polarized in y
     % - In the center of a cell if source polarized in z
     %
     % The switch(pola_source) takes care of the specific adjustment 
     % for avery polarisation
     
     sxGrid = round(sx./delta_x_max).*delta_x_max;
     szGrid = round(sz./delta_z_max).*delta_z_max;
     switch(pola_source)
         case 1
        % Source polarized in Z direction. The sources are located
        % on the corner of the grid cells.
        % We use the biggest grid_size to set the domain and the 
        % position of the sources. (see sxGrid and szGrid just above  
        % the switch(pola_source)
        %
        % Even if the grid size is reduced, the sources stay in a
        % corner of a cell and the domain size is not changed. The
        % top left corner is always at x0 and z0.
             for n = 1:NbSources
                 sx_idex(:,n) = round((sxGrid(n)-x0)./delta_x);
                 sz_idex(:,n) = round((szGrid(n)-z0)./delta_z);
             end
         case 0
        % Source polarized in Y direction. The sources are located
        % on the middle of the edge of a cell (in z direction).
        % We use the biggest grid_size to set the domain and the 
        % position of the sources. (see sxGrid and szGrid just above  
        % the switch(pola_source)
        
        % Every time the delta_z is reduced,the top left corner position z0
        % is adjusted in order that the sources stay at the same position
             for n = 2:NbGrid
                 if (~(delta_z(n)==delta_z(n-1)))
                     z0(n) = z0(n) + delta_z(n)/2;
                 else
                     z0(n) = z0(n-1);
                 end
             end
             for n = 1:NbSources % Adjust szGrid so sources ar located
                                 % in the midle of an edge in z direction
                 if ((szGrid(n)-sz(n))>0)
                     szGrid(n) = szGrid(n) - delta_z_max/2;
                 else
                     szGrid(n) = szGrid(n) + delta_z_max/2;
                 end
                 % calculate the idex for the position of sources 
                 sx_idex(:,n) = round((sxGrid(n)-x0)./delta_x);
                 sz_idex(:,n) = round((szGrid(n)-z0)./delta_z);
             end
         case 2
        % Source polarized in X direction. The sources are located
        % in the center of a cell.
        % We use the biggest grid_size to set the domain and the 
        % position of the sources. (see sxGrid and szGrid just above  
        % the switch(pola_source)
        
        % Every time the delta_z is reduced, the top left corner position 
        % x0 and z0 are adjusted in order that the sources stay at the 
        % same position
             for n = 2:NbGrid 
                 if (~(delta_z(n)==delta_z(n-1)))
                     z0(n) = z0(n) + delta_z(n)/2;
                     x0(n) = x0(n) + delta_x(n)/2;
                 else
                     z0(n) = z0(n-1);
                     x0(n) = x0(n-1);
                 end
             end
             for n = 1:NbSources % Adjust szGrid so sources ar located
                                 % in the center of a cell
                 if ((szGrid(n)-sz(n))>0)
                     szGrid(n) = szGrid(n) - delta_z_max/2;
                 else
                     szGrid(n) = szGrid(n) + delta_z_max/2;
                 end
                 if ((sxGrid(n)-sx(n))>0)
                     sxGrid(n) = sxGrid(n) - delta_x_max/2;
                 else
                     sxGrid(n) = sxGrid(n) + delta_x_max/2;
                 end
                % calculate the idex for the position of sources 
                sx_idex(:,n) = round((sxGrid(n)-x0)./delta_x);
                sz_idex(:,n) = round((szGrid(n)-z0)./delta_z); 
             end
     end  
end


function [epsilonMatrix_Grid, sigmaMatrix_Grid, muMatrix_Grid,Nx, Nz] = ...
                            Interpol(lx,lz,lx_Grid,lz_Grid,x0,z0,...
                            delta_x,delta_z,epsilonMatrix,muMatrix,...
                            sigmaMatrix,method)
    
    NbGrid = length(delta_x);
    temp = size(epsilonMatrix);
    NxPar = temp(2);
    NzPar = temp(1);
    Nx = round(lx_Grid./delta_x);
    Nz = round(lz_Grid./delta_z);
    
    for n = 1:NbGrid
        epsilonMatrix_Grid{n} = zeros(int16(Nz(n)),int16(Nx(n))); %#ok<*AGROW>
        sigmaMatrix_Grid{n} = zeros(int16(Nz(n)),int16(Nx(n)));
        muMatrix_Grid{n}(:,:) = zeros(int16(Nz(n)),int16(Nx(n)));
        %center of initial grid
        x = lx/NxPar/2:lx/NxPar:lx-lx/NxPar/2;
        z = lz/NzPar/2:lz/NzPar:lz-lz/NzPar/2;
        [X,Z] = meshgrid(x,z);
        % center of new grid
        x = (delta_x(n)/2:delta_x(n):lx_Grid-delta_x(n)/2) + x0(n);
        z = (delta_z(n)/2:delta_z(n):lz_Grid-delta_z(n)/2) + z0(n);
        [XGrid,ZGrid] = meshgrid(x,z);
        epsilonMatrix_Grid{n}(:,:) = ...
            interp2(X,Z,epsilonMatrix,XGrid,ZGrid,method);
        sigmaMatrix_Grid{n}(:,:) = ...
            interp2(X,Z,sigmaMatrix,XGrid,ZGrid,method);
        muMatrix_Grid{n}(:,:) = ...
            interp2(X,Z,muMatrix,XGrid,ZGrid,method);
    end

    %%% Now fill the nan's %%%%%%%
    for n = 1:NbGrid
        nanLocations = isnan(epsilonMatrix_Grid{n});
        nonNanLinearIndexes = find(~nanLocations);
        [zGood, xGood] = ...
            ind2sub([Nz(n),Nx(n)], nonNanLinearIndexes);
        xi = min(xGood);
        xf = max(xGood);
        zi = min(zGood);
        zf = max(zGood);
        %fill the corners
        epsilonMatrix_Grid{n}(1:zi,1:xi) =  ...
            epsilonMatrix_Grid{n}(zi,xi);
        epsilonMatrix_Grid{n}(zf:Nz(n),1:xi) =  ...
            epsilonMatrix_Grid{n}(zf,xi);
        epsilonMatrix_Grid{n}(1:zi,xf:Nx(n)) =  ...
            epsilonMatrix_Grid{n}(zi,xf);
        epsilonMatrix_Grid{n}(zf:Nz(n),xf:Nx(n)) = ...
            epsilonMatrix_Grid{n}(zf,xf);
        
        sigmaMatrix_Grid{n}(1:zi,1:xi) =  ...
            sigmaMatrix_Grid{n}(zi,xi);
        sigmaMatrix_Grid{n}(zf:Nz(n),1:xi) =  ...
            sigmaMatrix_Grid{n}(zf,xi);
        sigmaMatrix_Grid{n}(1:zi,xf:Nx(n)) =  ...
            sigmaMatrix_Grid{n}(zi,xf);
        sigmaMatrix_Grid{n}(zf:Nz(n),xf:Nx(n)) = ...
            sigmaMatrix_Grid{n}(zf,xf);
        
        muMatrix_Grid{n}(1:zi,1:xi) =  ...
            muMatrix_Grid{n}(zi,xi);
        muMatrix_Grid{n}(zf:Nz(n),1:xi) =  ...
            muMatrix_Grid{n}(zf,xi);
        muMatrix_Grid{n}(1:zi,xf:Nx(n)) =  ...
            muMatrix_Grid{n}(zi,xf);
        muMatrix_Grid{n}(zf:Nz(n),xf:Nx(n)) = ...
            muMatrix_Grid{n}(zf,xf);
         
         
        % fill the top and bottom
        for i = xi:xf
            epsilonMatrix_Grid{n}(1:zi,i) =  ...
                epsilonMatrix_Grid{n}(zi,i);
            epsilonMatrix_Grid{n}(zf:Nz(n),i) =  ...
                epsilonMatrix_Grid{n}(zf,i);
            
            sigmaMatrix_Grid{n}(1:zi,i) =  ...
                sigmaMatrix_Grid{n}(zi,i);
            sigmaMatrix_Grid{n}(zf:Nz(n),i) =  ...
                sigmaMatrix_Grid{n}(zf,i);
            muMatrix_Grid{n}(1:zi,i) =  ...
                muMatrix_Grid{n}(zi,i);
            muMatrix_Grid{n}(zf:Nz(n),i) =  ...
                muMatrix_Grid{n}(zf,i);
        end
         % fill the left and right
         for i = zi:zf
             epsilonMatrix_Grid{n}(i,1:xi) =  ...
                 epsilonMatrix_Grid{n}(i,xi);
             epsilonMatrix_Grid{n}(i,xf:Nx(n)) =  ...
                 epsilonMatrix_Grid{n}(i,xf);
             sigmaMatrix_Grid{n}(i,1:xi) =  ...
                 sigmaMatrix_Grid{n}(i,xi);
             sigmaMatrix_Grid{n}(i,xf:Nx(n)) =  ...
                 sigmaMatrix_Grid{n}(i,xf);
             muMatrix_Grid{n}(i,1:xi) =  ...
                 muMatrix_Grid{n}(i,xi);
             muMatrix_Grid{n}(i,xf:Nx(n)) =  ...
                 muMatrix_Grid{n}(i,xf);
         end
    end
         
end


function PlotModelAndGrid(lx,lz,epsilonMatrix,sigmaMatrix,sx,sz,rx,rz,...
                             lx_Grid,lz_Grid,Nx,Nz,epsilonMatrix_Grid,...
                             sigmaMatrix_Grid,sx_Grid,sz_Grid,delta_x,...
                             x0,z0,f,Npmlx,Npmlz)

    NbGrid = length(delta_x);
    Nf = length(f);

    temp = size(epsilonMatrix);
    NxPar = temp(2);
    NzPar = temp(1);

    idexIni = 1;
    for n = 1:NbGrid    
        f_min = f(idexIni,1);
        for i = idexIni:Nf
            if (f(i,2) > n)
                idexIni = i;
                f_max = f(idexIni,1);
                break
            end
            f_max = f(i,1);
        end
        X1 = lx/NxPar/2;
        X2 = lx-lx/NxPar/2;
        Z1 = lz/NzPar/2;
        Z2 = lz-lz/NzPar/2;
        
        figure(n)
        subplot(2,2,1)
        imagesc([X1 X2],[Z1 Z2],epsilonMatrix);
        xlabel('x')
        ylabel('z')
        title('\epsilon_r model');
        colorbar
        hold on;
        plot(sx,sz,'>k');
        for i = 1:length(sx)
            idex = int2str(i);
            text(sx(i)+lx/100,sz(i), strcat('S',idex),...
                                'Color','black','FontSize',14);
        end
        for i = 1:NxPar-1
            plot([i*lx/NxPar,i*lx/NxPar],[0,lz],'k-');
        end
        for j = 1:NzPar-1
            plot([0,lx],[j*lz/NzPar,j*lz/NzPar],'k-');
        end
       
    
        subplot(2,2,2)
        imagesc([X1 X2],[Z1 Z2],sigmaMatrix);
        xlabel('x')
        ylabel('z')
        title('\sigma model');
        colorbar
        hold on;
        plot(sx,sz,'>k');
        for i = 1:length(sx)
            idex = int2str(i);
            text(sx(i)+lx/100,sz(i), strcat('S',idex),...
                            'Color','black','FontSize',14);
        end
        for i = 1:NxPar-1
            plot([i*lx/NxPar,i*lx/NxPar],[0,lz],'k-');
        end
        for j = 1:NzPar-1
            plot([0,lx],[j*lz/NzPar,j*lz/NzPar],'k-');
        end
    
    
        % Next, we show the grid that is used to 
        % calculate the Green Function
     
        X1 = lx_Grid/Nx(n)/2;
        X2 = lx_Grid-lx_Grid/Nx(n)/2;
        Z1 = lz_Grid/Nz(n)/2;
        Z2 = lz_Grid-lz_Grid/Nz(n)/2;
        
        Z1 = Z1 + z0(n);
        Z2 = Z2 + z0(n);
        X1 = X1 + x0(n);
        X2 = X2 + x0(n);
        
        % To show the PML region, we copy the epsilon matrix
        % into a "temp" matrix and replace the values associated with
        % the PML regions by NaN
        temp  = epsilonMatrix_Grid{n}; 
        temp(:,:) = NaN;
        temp(Npmlz(n)+1:Nz(n)-Npmlz(n),Npmlx(n)+1:Nx(n)-Npmlx(n)) =  ...
            epsilonMatrix_Grid{n}(Npmlz(n)+1:Nz(n)-Npmlz(n),Npmlx(n)+1:Nx(n)-Npmlx(n));
        imAlpha=ones(size(temp));
        imAlpha(isnan(temp))=0;

        
        subplot(2,2,3);
        imagesc([X1 X2],[Z1 Z2],temp,'AlphaData',imAlpha);
        set(gca,'color',0.75*[1 1 1]);
        xlabel('x')
        ylabel('z')
        title({strcat('Grid size for \epsilon_r (Nx =',num2str(Nx(n)),';Nz=',num2str(Nz(n)),')');...
            strcat(num2str(f_min,2),'MHz < f <',num2str(f_max,2),'MHz')});
        hold on;
        plot(sx_Grid,sz_Grid,'*r');
        plot(sx,sz,'>k');
        plot(rx,rz,'ok','LineWidth',2);
        legend('Modified source position',...
                                    'Initial source position','Receptor');
        for i = 1:Nx(n)-1
            plot([i*lx_Grid/Nx(n)+x0(n),i*lx_Grid/Nx(n)+x0(n)],...
                                            [z0(n),lz_Grid+z0(n)],'k-');
        end
        for j = 1:Nz(n)-1
            plot([x0(n),lx_Grid+x0(n)],...
                [j*lz_Grid/Nz(n)+z0(n),j*lz_Grid/Nz(n)+z0(n)],'k-');
        end
       
        for i = 1:length(sx)
            idex = int2str(i);
            text(sx_Grid(i)+lx_Grid/100,sz_Grid(i), ...
                            strcat('S',idex),'Color','red','FontSize',14);
        end
         for i = 1:length(rx)
            idex = int2str(i);
            text(rx(i)+lx_Grid/50,rz(i), ...
                        strcat('R',idex),'Color','black','FontSize',14);
         end
        
        % To show the PML region, we copy the sigma matrix
        % into a "temp" matrix and replace the values associated with
        % the PML regions by NaN
        temp  = sigmaMatrix_Grid{n}; 
        temp(:,:) = NaN;
        temp(Npmlz(n)+1:Nz(n)-Npmlz(n),Npmlx(n)+1:Nx(n)-Npmlx(n)) = ...
                sigmaMatrix_Grid{n}(Npmlz(n)+1:Nz(n)-Npmlz(n),Npmlx(n)+1:Nx(n)-Npmlx(n));
        
        subplot(2,2,4);
        imagesc([X1 X2],[Z1 Z2],temp,'AlphaData',imAlpha);
        set(gca,'color',0.75*[1 1 1]);
        xlabel('x')
        ylabel('z')
        title({strcat('Grid size for \sigma (Nx =',num2str(Nx(n)),';Nz=',num2str(Nz(n)),')');strcat(num2str(f_min,2),...
                                    'MHz < f <',num2str(f_max,2),'MHz')});
       % colorbar
        hold on;
        plot(sx_Grid,sz_Grid,'*r');
        plot(sx,sz,'>k');
        plot(rx,rz,'ok','LineWidth',2);
        legend('Modified source position','Initial source position',...
                                                            'Receptor');
         for i = 1:Nx(n)-1
            plot([i*lx_Grid/Nx(n)+x0(n),i*lx_Grid/Nx(n)+x0(n)],...
                                               [z0(n),lz_Grid+z0(n)],'k-');
        end
        for j = 1:Nz(n)-1
            plot([x0(n),lx_Grid+x0(n)],[j*lz_Grid/Nz(n)+z0(n),...
                                              j*lz_Grid/Nz(n)+z0(n)],'k-');
        end
        for i = 1:length(sx)
            idex = int2str(i);
            text(sx_Grid(i)+lx_Grid/100,sz_Grid(i), strcat('S',idex),...
                                              'Color','red','FontSize',14);
        end
        for i = 1:length(rx)
            idex = int2str(i);
            text(rx(i)+lx_Grid/50,rz(i), strcat('R',idex),...
                                            'Color','black','FontSize',14);
        end 
    end
end