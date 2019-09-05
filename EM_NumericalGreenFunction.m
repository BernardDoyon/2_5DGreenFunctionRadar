% ======================================================================
% Copyright (c) June 2019, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================
function EM_NumericalGreenFunction(FileNameModelInfo,FileNameOut,...
                        FlagCoeffOpt)
                    
   % This function solves equation 6 and 7
   %
   % For every wavenumber Ky, the sparse matrix is constructed 
   % (with RadarSparseMatrixOpt_mex or RadarSparseMatrix_mex)
   % (see line 175, 188 (for Ky = 0) and 234, 246 (for the other values
   % of Ky)
   %
   % Then, equation 6 is solved with the backslash operation
   % (see line 197 (for Ky = 0) and 255 (for the other values
   % of Ky)
   %
   % Finally, the summation is calculated (line 260, we use odd and even 
   % property of Green's functions) until convergence criteria is
   % satisfied (see line 291)
   %    
   % The final Green function is stored for every cell of the grid
   %  
   %    Ex(Nf,1:Nz*Nx) (x component)
   %    Ez(Nf,1:Nz*Nx) (z component)
   %    Ey(Nf,1:Nz*Nx) (y component)
   %
   %    Example:
   %            Ex(22,3) is the complex Green function for the x componant 
   %            at frequency index "23" and at grid cell position index "3"
   %
   %    Note: The top left cell index is equal to 1. The value of the index
   %          is first increased from left to right and then from top to
   %          bottum
   %
  
    mu0 = 4*pi*1e-7;
    epsilon0 = 1/36/pi*1e-9;
     
    load(FileNameModelInfo, 'f','q','pola_source','pola_receptor',...
               'delta_x','delta_z','Nx', 'Nz','x0','z0','Npmlx','Npmlz',...
               'sx_idex', 'sz_idex','rx','rz','ry',...
               'epsilonMatrix_Grid','sigmaMatrix_Grid',...
               'muMatrix_Grid');
    
    
    omega_vect = 2*pi*f(:,1)*1e6;    %#ok<*NODEF> angular frequency (rad/s)
    omega_vect_q = omega_vect ...    % angular frequency with small 
                        + 1i*2*pi*q*1e6;    % imaginary component (rad/s)
    
    if(FlagCoeffOpt)
        beta(1,1) = 0.7523;
        beta(1,2) = 0.9221;
    end
   
    temp = size(f);
    Nf = temp(1);
    omega_max = 2*pi*f(Nf,1)*1e6;
    flag_gn = 0;
    
    % Evaluate smallest distance source-recpetor
    sx = sx_idex(1,:)*delta_x(1)+x0(1);
    sz = sz_idex(1,:)*delta_z(1)+z0(1);
    Rmin = zeros(length(sx),length(rx));
    for n = 1:length(sx)
        Rmin(n,:) = sqrt((sx(n)-rx).^2 + (sz(n)-rz).^2 + ry^2);
    end
    Rmin = min(min(Rmin));
    
    flagCoeffUsed = zeros(1,Nf);
    
    for w = 1:Nf
        if (f(w,3) == 1)
            gn = f(w,2); % grid number index for this frequency
            if ~(flag_gn == gn) %new grid size
                flag_gn = gn;
               
                epsilon = epsilonMatrix_Grid{gn}*epsilon0; %#ok<*USENS>
                mu = muMatrix_Grid{gn}*mu0;
                sigma = sigmaMatrix_Grid{gn};
                
                rx_idex = round((rx - x0(gn))/delta_x(gn));
                rz_idex = round((rz - z0(gn))/delta_z(gn));
                
                XIx = zeros(Nz(gn),Nx(gn),4);
                XIz = zeros(Nz(gn),Nx(gn),4);
                % XIx(z,x,1) -> halfprevious;
                % XIx(z,x,2) -> current;
                % XIx(z,x,3) -> halfnext;
                % XIx(z,x,4) -> next;
                % same notation for XIz
            end
            
            Ex{w} = zeros(1,Nx(gn)*Nz(gn)); %Allocate memory for Ex field
            Ey{w} = zeros(1,Nx(gn)*Nz(gn)); %Allocate memory for Ey field
            Ez{w} = zeros(1,Nx(gn)*Nz(gn)); %Allocate memory for Ez field
            
            omega = omega_vect_q(w);
            
            Y = -1i*omega*epsilon + sigma;  % admittivity
            Z = -1i*omega*mu;               % impedivity
            
            kR = min(min((abs(real(sqrt(-Y.*Z))))));
            kR = kR*Rmin;
            
            % phase velocity
            Vp = ((epsilon.*mu/2).*(sqrt(1+sigma.^2./epsilon.^2/...
                real(omega)^2)+1)).^(-1/2);
            VpMax = max(max(Vp)); % maximum value for phase velocity
            VpMin = min(min(Vp)); % minimum value for phase velocity
            
            Kc = real(omega)/VpMin; % When Ky > Kc => evanescent waves
            % for every cell
            
            % Set a limit to Ky
            % The summation should stop around Ky = Kc
            KyMax = 10*Kc;
            
            SlownessGhost=1/VpMax/1.2;
            dKy = omega_max/Nf*SlownessGhost/1.1;
            Nk = round(KyMax/dKy);
            Ky = KyMax*linspace(0,1,Nk+1); % ky de 0 à ky_max
            YspaceMax = 2*pi/(Ky(2)-Ky(1));
            
            S = zeros(3*Nx(gn)*Nz(gn),1); % allocate memory for source vector
            
            if(FlagCoeffOpt && kR > 10)
                fprintf('\n \nSolving for f = %5.2f MHz with optimal coefficients \n', f(w,1));
                flagCoeffUsed(w) = 1;
            else
                fprintf('\n \nSolving for f = %5.2f MHz with standard coefficients \n', f(w,1));
                flagCoeffUsed(w) = 0;
            end
            
            ETilde = zeros(Nk,3*Nx(gn)*Nz(gn));     % E field in k space
            % (all components)
            ETilde_x = zeros(Nk,Nx(gn)*Nz(gn));     % Ex in k space
            ETilde_y = zeros(Nk,Nx(gn)*Nz(gn));     % Ey in k space
            ETilde_z = zeros(Nk,Nx(gn)*Nz(gn));     % Ez in k space
            
            Ey_idex = 3*linspace(0,Nx(gn)*Nz(gn),Nx(gn)*Nz(gn)+1);
            Ey_idex = Ey_idex(2:Nx(gn)*Nz(gn)+1);  % Ey tilde index
            % Note:
            %   Ex tilde index = Ey_idex -2
            %   Ez tilde index = Ey_idex -1
            % set the index of source vector where the source is not zero
            idex_source = 3*(sz_idex(gn,:)-1)*Nx(gn) + ...
                3*sx_idex(gn,:) - pola_source;
            
            S(idex_source,1) = -1/delta_x(gn)/delta_z(gn);
            
            % Evaluate solution for Ky = 0
            % (we use even and odd property of Green function)
            %  Ky = 0 is the first term of summation
            % (see Ellefsen09) 
           
            % calculating coefficients for pml cells
            p = 3.5;
            for jj = 1:Nx(gn)
                for kk = 1:Nz(gn)
                    K2_5D = sqrt(-Y(kk,jj)*Z(kk,jj));
                    Kr = real(K2_5D);
                    Ki = imag(K2_5D);
                    [XIx(kk,jj,1), XIx(kk,jj,2), XIx(kk,jj,3),...
                        XIx(kk,jj,4)] = xsi_coef(jj,Npmlx(gn),...
                        Nx(gn),delta_x(gn),Kr,Ki,p);
                    [XIz(kk,jj,1), XIz(kk,jj,2), XIz(kk,jj,3),...
                        XIz(kk,jj,4)] = xsi_coef(kk,Npmlz(gn),...
                        Nz(gn),delta_z(gn),Kr,Ki,p);
                end
            end
            % end of calculation of coefficients for pml cells
            if (FlagCoeffOpt && kR > 10) % can we use optimal coefficients?
                if exist('RadarSparseMatrixOpt_mex','file') == 3
                    [lin,col,ele] = RadarSparseMatrixOpt_mex(int32(Nx(gn)), ...
                        int32(Nz(gn)),Ky(1),...
                        delta_x(gn),delta_z(gn),Y,Z,XIx,...
                        XIz,beta);
                else 
                    [lin,col,ele] = RadarSparseMatrixOpt(int32(Nx(gn)), ...
                        int32(Nz(gn)),Ky(1),...
                        delta_x(gn),delta_z(gn),Y,Z,XIx,...
                        XIz,beta);
                end
            else
                if exist('RadarSparseMatrix_mex','file') == 3
                    [lin,col,ele] = RadarSparseMatrix_mex(int32(Nx(gn)),...
                        int32(Nz(gn)),Ky(1),delta_x(gn),delta_z(gn),Y,Z,XIx,XIz);
                else
                    [lin,col,ele] = RadarSparseMatrix(int32(Nx(gn)),...
                        int32(Nz(gn)),Ky(1),delta_x(gn),delta_z(gn),Y,Z,XIx,XIz);                    
                end
            end
            
            A = sparse(lin,col,ele);
            ETilde(1,:)=A\S;
            if(pola_source) % z or x polarized
                Ex{w} = ETilde(1,Ey_idex-2)/YspaceMax; %#ok<*AGROW>
                Ez{w} = ETilde(1,Ey_idex-1)/YspaceMax;
                Ey{w} = 0;
            else % y polarized
                Ex{w} = 0; %#ok<*AGROW>
                Ez{w} = 0;
                Ey{w} = ETilde(1,Ey_idex)/YspaceMax;    
            end
            
            MyText = ['Approximative number k values for the sommation: '...
                num2str(round(Kc/dKy)) '. Solving for k = '];
            fprintf(MyText);
            
            for k = 2:Nk
                % calculating coefficients for pml cells
                p = 3.5;
                for jj = 1:Nx
                    for kk = 1:Nz
                        K2_5D = sqrt(-Y(kk,jj)*Z(kk,jj) - Ky(k)^2);
                        Ki = imag(K2_5D);
                        if (Ki<0)
                            Ki = -Ki;
                        end
                        Kr = real(K2_5D);
                        [XIx(kk,jj,1), XIx(kk,jj,2), ...
                            XIx(kk,jj,3), XIx(kk,jj,4)] = xsi_coef(jj,...
                            Npmlx(gn),Nx(gn),delta_x(gn),Kr,Ki,p);
                        [XIz(kk,jj,1), XIz(kk,jj,2),...
                            XIz(kk,jj,3), XIz(kk,jj,4)] = xsi_coef(kk,...
                            Npmlz(gn),Nz(gn),delta_z(gn),Kr,Ki,p);
                    end
                end
                % end of calculation of coefficients for pml cells
                if (FlagCoeffOpt && kR > 10)
                    if exist('RadarSparseMatrixOpt_mex','file') == 3
                        [lin,col,ele] = RadarSparseMatrixOpt_mex(int32(Nx(gn)), ...
                            int32(Nz(gn)),Ky(k),...
                            delta_x(gn),delta_z(gn),Y,Z,XIx,...
                            XIz,beta);
                    else
                        [lin,col,ele] = RadarSparseMatrixOpt(int32(Nx(gn)), ...
                            int32(Nz(gn)),Ky(k),...
                            delta_x(gn),delta_z(gn),Y,Z,XIx,...
                            XIz,beta);
                    end
                else
                    if exist('RadarSparseMatrix_mex','file') == 3
                        [lin,col,ele] = RadarSparseMatrix_mex(int32(Nx(gn)),...
                            int32(Nz(gn)),Ky(k),delta_x(gn),delta_z(gn),Y,Z,XIx,XIz);
                    else
                        [lin,col,ele] = RadarSparseMatrix(int32(Nx(gn)),...
                            int32(Nz(gn)),Ky(k),delta_x(gn),delta_z(gn),Y,Z,XIx,XIz);
                    end
                end
                A = sparse(lin,col,ele);
                
                ETilde(k,:)=A\S;
                ETilde_x(k,:) = ETilde(k,Ey_idex-2);
                ETilde_z(k,:) = ETilde(k,Ey_idex-1);
                ETilde_y(k,:) = ETilde(k,Ey_idex);
                
                if (pola_source) %source en x ou z
                    Ex{w} = Ex{w} + 2*ETilde_x(k,:)*cos(Ky(k)*ry)/YspaceMax;
                    Ez{w} = Ez{w} + 2*ETilde_z(k,:)*cos(Ky(k)*ry)/YspaceMax;
                    Ey{w} = Ey{w} + 2i*ETilde_y(k,:)*sin(Ky(k)*ry)/YspaceMax;
                else
                    Ex{w} = Ex{w} + 2i*ETilde_x(k,:)*sin(Ky(k)*ry)/YspaceMax;
                    Ez{w} = Ez{w} + 2i*ETilde_z(k,:)*sin(Ky(k)*ry)/YspaceMax;
                    Ey{w} = Ey{w} + 2*ETilde_y(k,:)*cos(Ky(k)*ry)/YspaceMax;
                end
                
                idex = linspace(0, (Nz(gn)-1)*Nx(gn), Nz(gn));
                idex = idex(rz_idex) + rx_idex;
                % idex = index of the cells where the receptors are
                switch pola_receptor
                    case 0
                        Gk_G = max(abs(ETilde_y(k,idex))./abs(Ey{w}(1,idex)));
                    case 1
                        Gk_G = max(abs(ETilde_z(k,idex))./abs(Ez{w}(1,idex)));
                    case 2
                        Gk_G = max(abs(ETilde_x(k,idex))./abs(Ex{w}(1,idex)));
                    case 3
                        Gk_G = max([max(abs(ETilde_y(k,idex))./abs(Ey{w}(1,idex))), ...
                            max(abs(ETilde_z(k,idex))./abs(Ez{w}(1,idex))), ...
                            max(abs(ETilde_x(k,idex))./abs(Ex{w}(1,idex)))]);
                end
                % Gk_G is the biggest relative correction to numerical solution
                % for all receptors. We only check the correction for the
                % polarisation of the receptor. If case 3, we check all
                % components and take the biggest correction to be Gk_G
                
                % we stop if the biggest relative correction is < 0.01
                if ((Ky(k) > Kc)&&(Gk_G<0.01))
                    break;
                end
                
                % show progression of the iterations on the command window
                if k>2
                    for j=1:length(TextToDisplay)
                        fprintf('\b'); % delete previous counter display
                    end
                end
                TextToDisplay = [num2str(k) ' (Gk/G = ' num2str(Gk_G,3) ')'];
                fprintf(TextToDisplay);
            end
        end
    end
       
    fprintf('\n');
    save(FileNameOut,'Ex','Ez','Ey','Nx','Nz','x0','z0','delta_x','delta_z','f',...
                        'q','sx_idex','sz_idex','pola_source',...
                        'pola_receptor','rx','rz','ry','flagCoeffUsed');
    
end
  


 
function [halfprevious, current, halfnext, next] = xsi_coef(idex, Npml,N,delta,Kr,Ki,p)

thick_left = Npml-1;
thick_right = Npml-0.5;

if (idex <= Npml) %left
    tmp = Npml-idex+1;
    dist_halfprev_l = tmp + 0.5;
    dist_current_l = tmp;
    dist_halfnext_l = tmp-0.5;
    dist_next_l =  tmp-1;
    
    T = thick_left*delta;
    Am = -(p+1)*(log(10^-8) + 2*T*Ki)/2/T/(Ki + Kr^2/Ki);
    
    halfprevious = 1 + (1+1i*Kr/Ki)*Am*(dist_halfprev_l/thick_left)^p;
    current = 1 + (1+1i*Kr/Ki)*Am*(dist_current_l/thick_left)^p;
    halfnext = 1 + (1+1i*Kr/Ki)*Am*(dist_halfnext_l/thick_left)^p;
    next = 1 + (1+1i*Kr/Ki)*Am*(dist_next_l/thick_left)^p;
elseif (idex >= N-Npml+1) %right
    tmp = idex-(N-Npml) - 1;
    dist_halfprev_r = tmp ;
    dist_current_r = tmp + 0.5;
    dist_halfnext_r = tmp + 1;
    dist_next_r =  tmp + 1.5;
    
    T = thick_right*delta;
    Am = -(p+1)*(log(10^-8) + 2*T*Ki)/2/T/(Ki + Kr^2/Ki);
    
    halfprevious = 1 + (1+1i*Kr/Ki)*Am*(dist_halfprev_r/thick_right)^p;
    current = 1 + (1+1i*Kr/Ki)*Am*(dist_current_r/thick_right)^p;
    halfnext = 1 + (1+1i*Kr/Ki)*Am*(dist_halfnext_r/thick_right)^p;
    next = 1 + (1+1i*Kr/Ki)*Am*(dist_next_r/thick_right)^p;
else
    halfprevious = 1;
    current = 1;
    halfnext = 1;
    next = 1;
end
end



