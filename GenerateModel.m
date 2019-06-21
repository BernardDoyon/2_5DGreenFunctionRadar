function GenerateModel(flag,FileName)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple funcition to generate a parameter model for
% sigma, epsilonr and mu_r
%
% 
% flag = 1 => Homogeneous Model
% flag = 2 => Heterogeneous Model (Clay-Sand)
% flag = 3 => Gradiant Model
%
% Numerical values for the parameters can be changed in this first 
% switch loop
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    currentFolder = pwd;
    DirDataOut = strcat(currentFolder,'\Model\');
    mkdir(DirDataOut);

    switch flag
        case 1 % Homogeneous model
            epsilon_r = 9.0;
            mu_r = 1.0;
            sigma = 1.0e-3;
            Nx = 10;
            Nz = 10;
        case 2 % 3 layer model (Top-Middle-Bottom)
                 % First value in vector is for Top
                 % In this example:
                 %      epsilon_r = 20 (Clay 1) 
                 %      epsilon_r = 40 (Sand)
                 %      epsilon_r = 10 (Clay2)
                 % We set 
                 %      Nz = [N1 N2 N3]
                 %  => First N1 cells in z is for Sand
                 %  => the N2 next cells in z are for Clay
                 %  => and the last N3 cells in z is for Sand
            epsilon_r = [40 20 40];
            mu_r = [1.0 1.0 1.0];
            sigma = [0.5 0.001 0.5];
            Nz = [1 10 1]; 
            Nx = 2;
        case 3
            % Gradient model
            epsilon_r = [20 40];
            mu_r = [1 1];
            sigma = [0.01 0.1];
            Nx = 10;
            Nz = 8;
    end
    
%   In this second switch loop, we set values in the matrices that
%   hold the epsilon_r,mu_r and sigma parameters. 
%
%   epsilonMatrix: epsilon_r Matrix
%   sigmaMatrix: sigma Matrix
%   muMatrix: mu_r Matrix 

    switch flag
        case 1
            epsilonMatrix = zeros(Nz,Nx);
            muMatrix = zeros(Nz,Nx);
            sigmaMatrix = zeros(Nz,Nx);
            epsilonMatrix(:,:) = epsilon_r; %scalar for homogeneous model
            muMatrix(:,:) = mu_r;
            sigmaMatrix(:,:) = sigma; 
        case 2
            NzTot = sum(Nz);
            epsilonMatrix = zeros(NzTot,Nx);
            muMatrix = zeros(NzTot,Nx);
            sigmaMatrix = zeros(NzTot,Nx);
            epsilonMatrix(:,:) = epsilon_r(2);
            epsilonMatrix(1:Nz(1),:) = epsilon_r(1);
            epsilonMatrix(NzTot-Nz(3)+1:NzTot,:) = epsilon_r(3);
            muMatrix(:,:) = mu_r(2);
            muMatrix(1:Nz(1),:) = mu_r(1);
            muMatrix(NzTot-Nz(3)+1:NzTot,:) = mu_r(3);
            sigmaMatrix(:,:) = sigma(2);
            sigmaMatrix(1:Nz(1),:) = sigma(1);
            sigmaMatrix(NzTot-Nz(3)+1:NzTot,:) = sigma(3);
        case 3
            e_grad_x = (epsilon_r(2) - epsilon_r(1))/Nx;
            e_grad_z = (epsilon_r(2) - epsilon_r(1))/Nz;
            s_grad_x = (sigma(2) - sigma(1))/Nx;
            s_grad_z = (sigma(2) - sigma(1))/Nz;
            m_grad_x = (mu_r(2)- mu_r(1))/Nx;
            m_grad_z = (mu_r(2)- mu_r(1))/Nz;
            epsilonMatrix = zeros(Nz,Nx);
            sigmaMatrix = zeros(Nz,Nx);
            muMatrix = zeros(Nz,Nx);
            for i = 1:Nx
                for j = 1:Nz
                    epsilonMatrix(j,i) = epsilon_r(1) + (j-1)*e_grad_z + (i-1)*e_grad_x;
                    sigmaMatrix(j,i) = sigma(1) + (j-1)*s_grad_z + (i-1)*s_grad_x;
                    muMatrix(j,i) = mu_r(1) + (j-1)*m_grad_z + (i-1)*m_grad_x;
                end
            end 
    end
    save(strcat(DirDataOut,FileName),'epsilonMatrix','muMatrix','sigmaMatrix')
end


