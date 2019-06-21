% ======================================================================
% Copyright (c) June 2019, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================
function main()
%
%   This is the main file to calculate the Green function for radar sources
%
%   Before using this main file, you must build the MEX file associated
%   with the two following .m files:
%           1) RadarSparseMatrix.m  -> RadarSparseMatrix_mex
%           2) RadarSparseMatrixOpt.m -> RadarSparseMatrixOpt_mex
%
%   If the Mex file can't be build, it is still possible to use the main
%   file BUT you must edit the file EM_NumericalGreenFunction.m and change 
%   every call of RadarSparseMatrix_mex by RadarSparseMatrix and every call
%   of RadarSparseMatrixOpt_mex by RadarSparseMatrixOpt. The programm will 
%   be much slower without the MEX file. 
%     
%   From an inital parameter model, the main file performs 3 operations.
%   The parameter model file must be in the "CurrentFolder\Model" directory 
%   and must contain 3 matrices for each of the parameters: 
%   epsilon_r (relative permittivity), mu_r (relative permeability)
%   and sigma (conductivity) on the domain. From this parameter model file 
%   (with the generic name "FileNameModel.mat", the main file can perform 
%   3 operations:
%
%   I) It will create a FileNameModelForwInfo.mat file associated 
%      with the forward modeling to perform. This file is created with 
%      the function InitializeForwardModeling(...) and will be located in
%      the CurrentFolder\Model directory. 
%   
%   II) It will call the function EM_NumericalGreenFunction(...)
%       to calculate the Green function everywhere on the numerical grid. 
%       All the info for this calculation are stored in the 
%       FileNameModelForwInfo.mat file created in step I) above.
%       The solution will be store in the FileNameGrennFunction file
%       located in CurrentFolder\Results. Interpolation is then used to 
%       calculate the solution at the receptors position. The solution at
%       receptors position will be store in FileNameGrennFunction_R file in
%       the variables:
%           ExR
%       
%
%   III) For the homgeneous model and the Three-Layers model, analytical 
%       solution for the Green function can be calculated and compare with 
%       the numerical solution
%
%
%   To create the FileNameModelForwInfo.mat file, the user has to 
%   provide informations about the simulation to perform. The informations 
%   are:
%
%   1) FileNameModel: 
%
%       FileNameModel is a .mat file with the 3 matrices
%       defining the model. The 3 matrices are:
%
%       1) epsilonMatrix(NzCell,NxCell):  relative permittivity
%       2) muMatrix(NzCell,NxCell)        relative permeability
%       3) sigmaMatrix(NzCell,NxCell)     conductivity  (Siemens/meter)
%
%       with
%           NzCell: Number of cell in the z direction
%           NxCell: Number of cell in the x direction
%       
%
%      NOTE: The FileNameModel.mat file has to be stored in the  
%            CurrentFolder\Model directory
%
%       All the other parameters of the model (for instance the length
%       of the domain or the position of the sources) are specified
%       in this main file. These parameters are: 
%
%   2) Nf:  Number of positive frequencies obtained from the 
%           Fourier transform of the intial time radar trace
%   3) fmax: Frequency max: Maximum value of the frequency (MHz) 
%   4) fmin: Frequency min  Minimum value of the frequency (MHz)
%
%   5) idexFrequency 
%                   Index of the frequencies for which a Green function has
%                   to be calculated.
%                   Example:
%                   Nf = 10; (10 positive frequencies) fmax = 100; fmin = 0
%                   The positive frequencies after the Fourier Transform of 
%                   the initial Radar trace are:
%                   f =[10 20 30 40 50 60 70 80 90 100] MHz
%                   If we set idexFrequency = [3 6 9], the Green function
%                   will be calculated for f = 30, 60 and 90 MHz. 
%        
%   6) q:    Imaginary component of frequency (MHz)  
%        Note: Rule of thumb:
%                (fmax-fmin)/2/Nf < q <(fmax-fmin)/Nf
%
%   7) lx:      x length of domain (meters)
%   8) lz:      z length of domain (meters)
%   9) sx:      x position of sources (meters)
%   10) sz:      z position of sources (meters)
%   11) rx:     x position of receptors (meters)
%   12) rz:     z position of receptors (meters)
%   13) ry:     y position of receptors (meters)
%   14) lpmlx:  x length of pml region (meters)
%   15) lpmly:  y length of pml region (meters)
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
%           of the  
%           - "y" component of the Green function (pola_receptor = 0), 
%           - "z" component of the Green function (pola_receptor = 1), 
%           - "x" component of the Green function (pola_receptor = 2), 
%           - all components of the Green function (pola_receptor = 3).
%  
%   18) delta_x: Grid spacing in the x direction (meters)
%
%   19) delta_z: Grid spacing in the z dimension (meters)
%
%       Note:
%           If delta_x and delta_z = 0, a grid spacing will be calculated 
%           in order that the minimum number of points per wavelength 
%           requested by the user is respected 
%           (see variable "NbGridPointPerWavelength" below). 
%           For delta_x and delta_z = 0, grid spacing can be different from 
%           one frequency to the other. 
%           A value of delta_x and delta_z are first calculated for the 
%           highest frequency. 
%           (to respect the minimum number of points per wavelength AND to 
%           be a divider of the initial parameter model grid size) 
%           The value of delta_x and delta_z are double every 
%           time the criteria for the NbGridPointPerWavelength is 
%           respected. It is possible to restrict the number of different 
%           Delta values (see variable NbGridSizes below)
%
%           The value for delta_x and delta_z CAN'T be smaller than the
%           initial parameter model grid size.
%
% 20) NbGridPointPerWavelength: The minimum number of points per wavelength
%                               to be respected for a given frequency.
%   
% 21) NbGridSizes: The maximum nuber of different Delta values. 
%
% 22) InterpolationMethod:  The interpolation method for the parameter 
%                           model. Choices are:
%                                  - 'linear',
%                                  - 'nearest',
%                                  - 'cubic',
%                                  - 'spline'
%       Note:   
%           The numerical grid to calculate the Green's functions will 
%           usually be smaller than the inital parameter model.  The 
%           parameter epsilon_r, mu_r and sigma need to be 
%           interpolated on this smaller new grid. The   
%
% 23)*flagPlotModel: When this flag is set to 1, a color image of the
%                   initial model for epsilon and sigma are plot for the
%                   initial domain. Color images of the interpolated
%                   models for epsilon and sigma are also plot for each
%                   numerical grid.   
%
%  
% To calculate the Green function, we can use the standar coefficients or
% the optimal coefficients. We use a flag to indicate the coefficients 
% 
%       flagCoeff: % 1 = Optimal Coefficient; 0 = Standard coefficients
%
%   Note:   The optimal coefficients are no use for kr < 10 (approx) where
%           k is the wavenumber and r is the distance from the source. For 
%           each frequency, we check the ccondition kr < 10. The standard 
%           coefficent will be used if the condition is not satisfied. 
%
%   Some models have been preset with a "flagModel" parameter
%   (see the switch(flagModel) below):
%
%           flagModel = 1 -> Homogeneous Model 
%           flagModel = 2 -> three-layers model
%           flagModel = 3 -> Gradient Model
%           flagModel = 4 -> Your Model (TO COMPLETE!)
%
%    Note: To modify the values for epsilon_r, mu_r and sigma in the 
%           Homogeneous, three-layers or gradient model, edit the 
%           GenerateModel.m file a change the values for these parameters 
%           in the file.
%


currentFolder = pwd;
DirDataOut = strcat(currentFolder,'\Results\');
mkdir(DirDataOut);

flagCoeff = 1; % 1 = Optimal Coefficient; 0 = Standard coefficients

flagModel = 1;      % 1-HomogeneousModel; 2-Three-Layer Medium
                    % 3-GradientModel 4-Your Model
                    
switch(flagModel)
    case 1
        FileNameGrennFunction = 'GreenFunctionsHomogeneous';
        FileNameModel = 'HomogeneousModel';
        flagPlotModel = 1;
        Nf = 45;                % Number of positive frequencies obtained 
                                % from the Fourier transform of the intial 
                                % time radar trace
        fmax = 150;             % Maximum value of the frequency 
        fmin = 0;               % Minimum value of the frequency
        idexFrequency = linspace(1,15,15)*3;
                                % The Green function will be calculated
                                % only for the frequency specified by 
                                % idexFrequency         
        q = 5;
        lx = 4;                 % x length of domain in meters
        lz = 1;                 % z length of domain in meters
        sx = [0.5 0.5];            % x position of sources in meters
        sz = [0.5 0.7];            % z position of sources in meters
        rx = [3.52 3.52 3.52];            % x position of receptors in meters
        rz = [0.3 0.55 0.7]; % z position of receptors in meters
        ry = -0.1;              % y position of receptors in meters
        lpmlx = 0.4;
        lpmlz = 0.4;
        pola_source = 1;
        pola_receptor = 3;
        delta_x = 0.04;% 0.1;
        delta_z = 0.04;%.1;
        NbGridPointPerWavelength = 20;
        NbGridSizes = 2;
        InterpolationMethod = 'nearest'; 
        GenerateModel(flagModel,FileNameModel);
    case 2
        FileNameGrennFunction = 'GreenFunctionsSandClay';
        FileNameModel = 'SandClayModel';
        flagPlotModel = 0;
        Nf = 24;
        fmax = 300;
        fmin = 0;
        idexFrequency = linspace(1,4,4)*4;
        q = 12.5;
        lx = 0.8;
        lz = 0.8;
        sx = [0.1 0.1];
        sz = [0.4 0.3];
        rx = [0.7 0.7 0.7];
        rz = [0.05 0.3 0.75];
        ry = -0.1;
        lpmlx = 0.1;
        lpmlz = 0.1;
        pola_source = 1;
        pola_receptor = 1;
        delta_x = 0;%0.01;
        delta_z = 0;%0.01;
        NbGridPointPerWavelength = 16;
        NbGridSizes = 3;
        InterpolationMethod = 'nearest'; 
        GenerateModel(flagModel,FileNameModel);
    case 3
        FileNameGrennFunction = 'GreenFunctionsGradient';
        FileNameModel = 'GradientModel';
        flagPlotModel = 1;
        Nf = 24;
        fmax = 150;
        fmin = 0;
        idexFrequency =  linspace(1,6,6)*4;
        q = 5;
        lx = 4;                 
        lz = 2;                
        sx = [1 1];            
        sz = [1 1.1];            
        rx = [3 3 3 3 3 3 3];      
        rz = [0.4 0.6 0.8 1 1.2 1.4 1.6];    
        ry = -0.1;          
        lpmlx = 0.4;
        lpmlz = 0.4;
        pola_source = 2;
        pola_receptor = 1;
        delta_x = 0;
        delta_z = 0;
        NbGridPointPerWavelength = 10;
        NbGridSizes = 3;
        InterpolationMethod = 'nearest';
        GenerateModel(flagModel,FileNameModel);
    case 4
        return; 
       
        % uncomment when FileNameModel is created for your model 
        % and complete the values for all the parameters
        
        %{
        FileNameGrennFunction = 'YourModelGreenFunction';
        FileNameModel = 'YourModel';
        Nf = ;
        fmax = ;
        fmin = ;
        idexFrequency;
        q = ;
        lx = ;
        lz = ;
        sx = ;
        sz = ;
        rx = ;
        rz = ;
        ry = ;
        lpmlx = ;
        lpmlz = ;
        pola_source = ;
        pola_receptor = ;
        delta_x = ;
        delta_z = ;
        NbGridPointPerWavelength = ;
        NbGridSizes = ;
        InterpolationMethod = ;
        %}
end

FullFileNameModel = strcat(currentFolder,'\Model\',FileNameModel);
FullFileNameModelInfo = strcat(FullFileNameModel,'InfoForwMod');

if (flagCoeff == 1) % Use optimal coefficients if appropriate
    FileNameNum = strcat(FileNameGrennFunction,'Opt');
else
    FileNameNum = strcat(FileNameGrennFunction,'Std');
end


InitializeForwardModeling(FullFileNameModel,Nf,fmin,fmax,idexFrequency,q,lx,lz,...
                           sx,sz,rx,rz,ry,lpmlx,lpmlz,pola_source,...
                           pola_receptor,delta_x,delta_z,NbGridPointPerWavelength,...
                           NbGridSizes,InterpolationMethod,flagPlotModel);

if (flagPlotModel)   
    message = 'Enlarge the figure to see if it corresponds to the';
    message = strcat(message,' model and then click ok');
    waitfor(msgbox(message));
    close;
end

answer = questdlg('Calculate Green function?','Yes','No');

if (strcmp(answer,'Yes'))
%%% Calculate the Green function everywhere on the grid %%%%%%%%%%%%%
FullFileNameNum = strcat(DirDataOut,FileNameNum);
EM_NumericalGreenFunction(FullFileNameModelInfo,...
                            FullFileNameNum,...
                              flagCoeff);

%%%%%%%%  Interpolation of numerical solution at receptors position %%%%%
load(FullFileNameNum,'Ex','Ez','Ey','Nx','Nz','x0','z0',...
                                        'delta_x','delta_z','f','sx_idex','sz_idex',...
                                        'pola_source','rx','rz','ry',...
                                        'flagCoeffUsed');                                    
[ExR, EzR, EyR] = Interpol(f,Nx,Nz,x0,z0,delta_x,delta_z,Ex,Ez,Ey,rx,rz); %#ok<*ASGLU>
save(strcat(DirDataOut,FileNameNum,'_R'),'f','ExR','EzR','EyR','delta_x',...
                                 'delta_z','sx_idex','sz_idex','rx','rz',...
                                 'ry','pola_source','flagCoeffUsed');
end


if (flagModel < 3)
    answer = questdlg('Compare with analytical expression?','Yes','No');

    if (strcmp(answer,'Yes'))
    %%%%%%%  Calculate analytical solution if available %%%%%%%%%%%%%%%%%%
        if (flagModel == 1) 
        % Calculate analytical Green function if model is
        % homogeneous and compare with numerical solution
        % load model info
        load(FullFileNameModelInfo, 'ry','rx','rz','delta_z','delta_x',...
            'sx_idex','sz_idex','f','q',...
            'Nx','Nz','x0','z0','pola_source',...
            'epsilonMatrix_Grid',...
            'sigmaMatrix_Grid',...
            'muMatrix_Grid');
        
        FileNameAnalytic = strcat(FileNameGrennFunction,'Ana_R');
        FileNameNum0 = strcat(FileNameGrennFunction,'Std_R');
        FileNameNum1 = strcat(FileNameGrennFunction,'Opt_R');
        
        % calculate position of sources
        switch(pola_source) 
            case 0 % source polarized in y
                 sx = sx_idex(1,:)*delta_x(1)+x0(1);  %#ok<*NODEF>
                 sz = sz_idex(1,:)*delta_z(1)+z0(1)-delta_z(1)/2;
            case 1 % source polarized in z
                sx = sx_idex(1,:)*delta_x(1)+x0(1);  %#ok<*NODEF>
                sz = sz_idex(1,:)*delta_z(1)+z0(1);
            case 2 % source polarized in x
                sx = sx_idex(1,:)*delta_x(1)+x0(1)-delta_x(1)/2;  %#ok<*NODEF>
                sz = sz_idex(1,:)*delta_z(1)+z0(1)-delta_z(1)/2;
        end
        
        % get values for epsilon_r, mu_r and sigma
        % for homogeneous model
        mu_r = muMatrix_Grid{1}(1,1);  %#ok<*USENS>
        epsilon_r = epsilonMatrix_Grid{1}(1,1);
        sigma = sigmaMatrix_Grid{1}(1,1);
        FullFileNameAnalytic = strcat(DirDataOut,FileNameAnalytic);
        AnalyticalSolutionReceptorsHomogeneous(mu_r,epsilon_r,sigma,...
                    sx,sz,f,q,rx,rz,ry,pola_source,FullFileNameAnalytic);
        switch (flagCoeff)
            case 0
                FullFileNameNum = strcat(DirDataOut,FileNameNum0);
            case 1
                FullFileNameNum = strcat(DirDataOut,FileNameNum1);
        end
        GreenFunctionErrorReceptors(FullFileNameAnalytic, FullFileNameNum);   
       % GreenFunctionErrorReceptorsComp(FullFileNameAnalytic, ...
       %                             strcat(DirDataOut,FileNameNum0),...
       %                             strcat(DirDataOut,FileNameNum1));
    elseif (flagModel == 2) % Calculate analytical Green function for
        % Sand-Clay model and compare with numerical
        % solution
        load(FullFileNameModelInfo, 'ry','rx','rz','delta_x','delta_z',...
            'sx_idex','sz_idex','f','q',...
            'Nx','Nz','x0','z0',...
            'epsilonMatrix_Grid',...
            'sigmaMatrix_Grid',...
            'muMatrix_Grid');
        
        FileNameAnalytic = strcat(FileNameGrennFunction,'Ana_R');
        FileNameNum0 = strcat(FileNameGrennFunction,'Std_R');
        FileNameNum1 = strcat(FileNameGrennFunction,'Opt_R');
        % calculate position of sources
        % (analytical solution coded for source polarized in z)
        sx = sx_idex(1,:)*delta_x(1) + x0(1);
        sz = sz_idex(1,:)*delta_z(1) + z0(1);
        
        % get values for epsilon_r, mu_r and sigma
        % for Sand-Clay model.
        epsilon_r = epsilonMatrix_Grid{1}(:,1);
        sigma = sigmaMatrix_Grid{1}(:,1);
        mu_r = muMatrix_Grid{1}(:,1);
        
        AnalyticalSolutionReceptorsThreeLayers(mu_r,epsilon_r, sigma,...
            delta_z(1),z0(1),sx,sz,f,q,rx,rz,ry,pola_source,...
            strcat(DirDataOut,FileNameAnalytic));
        
        FullFileNameAnalytic = strcat(DirDataOut,FileNameAnalytic);
        switch (flagCoeff)
            case 0
                FullFileNameNum = strcat(DirDataOut,FileNameNum0);
            case 1
                FullFileNameNum = strcat(DirDataOut,FileNameNum1);
        end
        GreenFunctionErrorReceptors(FullFileNameAnalytic, FullFileNameNum);
       % GreenFunctionErrorReceptorsComp(FullFileNameAnalytic, ...
       %                            strcat(DirDataOut,FileNameNum0),...
       %                             strcat(DirDataOut,FileNameNum1));
        end
    end 
end

end


function [ExR, EzR, EyR] = Interpol(f,Nx,Nz,x0,z0,delta_x,delta_z,Ex,Ez,Ey,rx,rz)
        
    Nf = length(f(:,1));
    ExR = zeros(Nf,length(rx));
    EzR = zeros(Nf,length(rx));
    EyR = zeros(Nf,length(rx));

    Nf = size(Ex);
    Nf = Nf(2);
    
    for w = 1:Nf
        if (f(w,3) == 1)        
            gn = f(w,2);
            
            ExGrid{w} = zeros(Nz(gn),Nx(gn)); %#ok<*AGROW>
            EzGrid{w} = zeros(Nz(gn),Nx(gn));
            EyGrid{w} = zeros(Nz(gn),Nx(gn));
            if (~(isempty(Ex{w})))
                
                for n=1:Nz(gn)
                    ExGrid{w}(n,:) = Ex{w}(Nx(gn)*(n-1)+1:Nx(gn)*n);
                    EzGrid{w}(n,:) = Ez{w}(Nx(gn)*(n-1)+1:Nx(gn)*n);
                    EyGrid{w}(n,:) = Ey{w}(Nx(gn)*(n-1)+1:Nx(gn)*n);
                end
                
                % for z polarization
                x = (delta_x(gn):delta_x(gn):Nx(gn)*delta_x(gn)) + x0(gn);
                z = (delta_z(gn):delta_z(gn):Nz(gn)*delta_z(gn)) + z0(gn);
                [XGrid,ZGrid] = meshgrid(x,z);
                
                EzR(w,:) = interp2(XGrid,ZGrid,EzGrid{w},rx,rz,'linear');
                
                % for y polarization
                ZGrid = ZGrid - delta_z(gn)/2;
                EyR(w,:) = interp2(XGrid,ZGrid,EyGrid{w},rx,rz,'linear');
                
                % for x polarization (ZGrid already OK)
                XGrid = XGrid - delta_x(gn)/2;
                ExR(w,:) = interp2(XGrid,ZGrid,ExGrid{w},rx,rz,'linear');      
            end
        end
    end
    
end

