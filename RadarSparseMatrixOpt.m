% ======================================================================
% Copyright (c) June 2019, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================


function [lin, col, ele] = RadarSparseMatrixOpt(Nx, Nz,Ky,delta_x,...
                                            delta_z,Y,Z,XIx,XIz,beta)
 
	% This file creates the elements of the sparse matrix associated
    % with the finite difference problem to solve in 2.5 D(equation 6) with
    % the optimal coefficients.
    % 
    % Input (the size of the different variables are also given for
    %           compiling the MEX file associated with this function)
    %       
    %   Nx  [size: int32(1 x 1)]: Total number of grid points in x 
    %   Nz  [size: int32(1 x 1)]: Total number of grid points in z 
    %   delta_x [size: double(1 x 1)]: grid size in x
    %   delta_z [size: double(1 x 1)]: grid size in z
    %   Ky [size: double(1 x 1)]: Value of the wave number
    %   Y [size: complex(double(:inf x :inf)) admittivity matrix
    %   Z (size: complex(double(:inf x :inf)) impedivity matrix
    %   XIx (size: complex(double(:inf x :inf x 4)) 
    %                   absorbing coefficents for boundary in x
    %   XIz (size: complex(double(:inf x :inf x 4)) 
    %                   absorbing coefficents for boundary in z
    %   beta [size: double(1 x 2)] : Values of the optimal coefficents

    Ne = (Nx-2)*(Nz-2)*27 + (Nx-4)*(Nz-4)*16 + ...  
             + (Nx + Nz-4)*9 + (Nx+Nz-4)*12 + (Nx+Nz-6)*12 + 24 + Nz + Nx;
                                    %Number of non-null element
    lin=zeros(1,Ne); % index of line element
    col=zeros(1,Ne);  % index of colomn element
    ele=zeros(1,Ne)*(1 + 1i);% value of element
    idex = 0;
       
       
    C11 = zeros(3,3)*(1+1i);
    C12 = zeros(3,3)*(1+1i); %#ok<*NASGU>
    C13 = zeros(3,3)*(1+1i);
    C21 = zeros(3,3)*(1+1i);
    C22 = zeros(3,3)*(1+1i);
    C23 = zeros(3,3)*(1+1i);
    C31 = zeros(3,3)*(1+1i);
    C32 = zeros(3,3)*(1+1i);
    C33 = zeros(3,3)*(1+1i);
    c_index=2;


    % Usefull constant fot Engquist-Majda condition at boundary

    rx = 1.0 /delta_x/delta_x;
	rz = 1.0 /delta_z/delta_z;
    rxHalf = rx / 2.0;
	rzHalf = rz / 2.0 ;
    ky2Half = Ky*Ky/2.0 ;
    inv_dxz = 1/sqrt(delta_x^2+delta_z^2);
    
    Zinv = 1./Z;
    
    
    j=1;
    for k =1:Nz
        n1 = 3*(j + Nx*(k-1))-2; % Ex
        idex = idex+1;
        lin(idex) = n1;
        col(idex) = n1;
        ele(idex) = 1;
    end
        
    k=Nz;
    for j =1:Nx
        n2 = 3*(j + Nx*(k-1))-1; % Ez
        idex = idex+1;
        lin(idex) = n2;
        col(idex) = n2;
        ele(idex) = 1;
    end
    
    
    
   %%% Fill the non-empty terms for the corners of grid %%%%%%%%%%%%%%
    
    % Top left corner
    % Champ Ex (j =2; k=1); Ez (j=1;k=1); Ey (j=1;k=1); 
    
    idex_j = [2 1 1];
    idex_k = [1 1 1];
    is = 1; % indice s pour offset colonne
    it = 1; % indice t pour offset colonne
    m = 3*(is + Nx*it); % offset de la colonne à remplir
   
    n1 = 3*(idex_j(1) + Nx*(idex_k(1)-1))-2; % Ex
    n2 = 3*(idex_j(2) + Nx*(idex_k(2)-1))-1; % Ez
    n3 = 3*(idex_j(3) + Nx*(idex_k(3)-1));   % Ey   
    vect_lin = [n1 n2 n3];
         
    for ix =1:3
        Kw = sqrt(-Y(idex_k(ix),idex_j(ix))*Z(idex_k(ix),idex_j(ix))-Ky*Ky);
        idex = idex+1;
        lin(idex) = vect_lin(ix);
        col(idex) = vect_lin(ix);
        ele(idex) = -inv_dxz + 1i*Kw/2;
        idex = idex+1;
        lin(idex) = vect_lin(ix);
        col(idex) = vect_lin(ix)+m;
        ele(idex) = inv_dxz + 1i*Kw/2;
    end
    
    
      % Bottom left
    % Champ Ex (j=2; k=Nz); Ez (j=1;k=Nz-1); Ey (j=1;k=Nz); 
    
    idex_j = [2 1 1];
    idex_k = [Nz Nz-1 Nz];
    is = 1; % indice s pour offset colonne
    it = -1; % indice t pour offset colonne
    m = 3*(is + Nx*it); % offset de la colonne à remplir
   
    n1 = 3*(idex_j(1) + Nx*(idex_k(1)-1))-2; % Ex
    n2 = 3*(idex_j(2) + Nx*(idex_k(2)-1))-1; % Ez
    n3 = 3*(idex_j(3) + Nx*(idex_k(3)-1));   % Ey   
    vect_lin = [n1 n2 n3];
         
    for ix =1:3
        Kw = sqrt(-Y(idex_k(ix),idex_j(ix))*Z(idex_k(ix),idex_j(ix))-Ky*Ky);
        idex = idex+1;
        lin(idex) = vect_lin(ix);
        col(idex) = vect_lin(ix);
        ele(idex) = -inv_dxz + 1i*Kw/2;
        idex = idex+1;
        lin(idex) = vect_lin(ix);
        col(idex) = vect_lin(ix)+m;
        ele(idex) = inv_dxz + 1i*Kw/2;
    end
   
    % Top right 
    % Champ Ex (j=Nx; k=1); Ez (j=Nx;k=1); Ey (j=Nx;k=1); 
    
    idex_j = [Nx Nx Nx];
    idex_k = [1 1 1];
    is = -1; % indice s pour offset colonne
    it = 1; % indice t pour offset colonne
    m = 3*(is + Nx*it); % offset de la colonne à remplir
   
    n1 = 3*(idex_j(1) + Nx*(idex_k(1)-1))-2; % Ex
    n2 = 3*(idex_j(2) + Nx*(idex_k(2)-1))-1; % Ez
    n3 = 3*(idex_j(3) + Nx*(idex_k(3)-1));   % Ey   
    vect_lin = [n1 n2 n3];
         
    for ix =1:3
        Kw = sqrt(-Y(idex_k(ix),idex_j(ix))*Z(idex_k(ix),idex_j(ix))-Ky*Ky);
        idex = idex+1;
        lin(idex) = vect_lin(ix);
        col(idex) = vect_lin(ix);
        ele(idex) = inv_dxz - 1i*Kw/2;
        idex = idex+1;
        lin(idex) = vect_lin(ix);
        col(idex) = vect_lin(ix)+m;
        ele(idex) = -inv_dxz - 1i*Kw/2;
    end
   
      % Bottom right
    % Champ Ex (j=Nx; k=Nz); Ez (j=Nx;k=Nz-1); Ey (j=Nx;k=Nz); 
    
    idex_j = [Nx Nx Nx];
    idex_k = [Nz Nz-1 Nz];
    is = -1; % indice s pour offset colonne
    it = -1; % indice t pour offset colonne
    m = 3*(is + Nx*it); % offset de la colonne à remplir
   
    n1 = 3*(idex_j(1) + Nx*(idex_k(1)-1))-2; % Ex
    n2 = 3*(idex_j(2) + Nx*(idex_k(2)-1))-1; % Ez
    n3 = 3*(idex_j(3) + Nx*(idex_k(3)-1));   % Ey   
    vect_lin = [n1 n2 n3];
         
    for ix =1:3
        Kw = sqrt(-Y(idex_k(ix),idex_j(ix))*Z(idex_k(ix),idex_j(ix))-Ky*Ky);
        idex = idex+1;
        lin(idex) = vect_lin(ix);
        col(idex) = vect_lin(ix);
        ele(idex) = inv_dxz - 1i*Kw/2;
        idex = idex+1;
        lin(idex) = vect_lin(ix);
        col(idex) = vect_lin(ix)+m;
        ele(idex) = -inv_dxz - 1i*Kw/2;
    end
   
    
    %%% Étant donné la grille décalée, il faut faire attention aux arêtes
    %%% Par exemple pour l'arête F4 (droite), on peut évaluer les équations pour Ex
    %%% à l'aide des fonctions f_C11, f_C12 et f_C13 qui demandent de
    %%% connaître le champ à gauche de l'arête (soit le point j=Nx-1
    %%% Pour les autres équations Ey et Ez, on utilise la condition
    %%% Engquist-Majda
    
   %%%%%%%%%%%%%%% FRONTIÈRE EN HAUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   
   k = Nz; % Jy est traité par Engquist-Majda et Jx comme un coin
   j = 2;
		
   Kw = sqrt(-Y(k,j)*Z(k,j)); % on veut grandeur de K
   tmp1 = -Kw*Kw + ky2Half + 2.0*1i*Kw/delta_z + rx;
   tmp2 = -Kw*Kw + ky2Half - 2.0*1i*Kw/delta_z + rx;
   
   n3 = 3*(j + Nx*(k-1));      % ligne assoicée à Jy, Jx est un coin
   is = [-1 -1 0 0 1 1];
   it = [-1 0 -1 0 -1 0];
   
   C33(c_index-1,c_index-1) = -rxHalf;
   C33(c_index-1,c_index)= -rxHalf;
   C33(c_index,c_index-1)= tmp1;
   
   C33(c_index,c_index) = tmp2;
   C33(c_index+1,c_index-1) = -rxHalf;
   C33(c_index+1,c_index) = -rxHalf;

   for ix = 1:length(is)
       m = 3*(is(ix) + Nx*it(ix));
       idex = idex+1;
       lin(idex) = n3;
       col(idex) = n3 + m;
       ele(idex) = C33(c_index+is(ix),c_index+it(ix));
   end

   for j = 3:Nx-1 % k = Nz
       Kw = sqrt(-Y(k,j)*Z(k,j)); % on veut grandeur de K
       tmp1 = -Kw*Kw + ky2Half + 2.0*1i*Kw/delta_z + rx;
       tmp2 = -Kw*Kw + ky2Half - 2.0*1i*Kw/delta_z + rx;
       
       n1 = 3*(j + Nx*(k-1))-2;    % ligne assoiée à x
       n3 = 3*(j + Nx*(k-1));      % ligne assoicée à y
       is = [-1 -1 0 0 1 1];
       it = [-1 0 -1 0 -1 0];
       
       %mêmes équations pour x et y
       
       C11(c_index-1,c_index-1) = -rxHalf;
       C11(c_index-1,c_index)= -rxHalf;
       C11(c_index,c_index-1)= tmp1;
       
       C11(c_index,c_index) = tmp2;
       C11(c_index+1,c_index-1) = -rxHalf;
       C11(c_index+1,c_index) = -rxHalf;
     
       vect_lin = [n1 n3];
       for ix = 1:length(is)
           for idex_line = 1:2
               n = vect_lin(idex_line);
               m = 3*(is(ix) + Nx*it(ix));
               idex = idex+1;
               lin(idex) = n;
               col(idex) = n + m;
               ele(idex) = C11(c_index+is(ix),c_index+it(ix));
           end
       end
       
   end
   
   %%%%%%%%%%%%%%%%% FIN FRONTIÈRE EN HAUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
  
    %%%%%%%%%%%%%%% FRONTIÈRE AU BAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   k = 1; % Jy est traité par Engquist-Majda et Jx comme un coin
   j = 2;
   	
   Kw = sqrt(-Y(k,j)*Z(k,j)); % on veut grandeur de K
   tmp1 = Kw*Kw - ky2Half + 2.0*1i*Kw/delta_z - rx;
   tmp2 = Kw*Kw - ky2Half - 2.0*1i*Kw/delta_z - rx;
   
   n3 = 3*(j + Nx*(k-1));      % ligne assoicée à Jy, Jx est un coin
   is = [-1 -1 0 0 1 1];
   it = [0 1 0 1 0 1];
   
   C33(c_index-1,c_index) = rxHalf;
   C33(c_index-1,c_index+1)= rxHalf;
   C33(c_index,c_index)= tmp1;
  
   C33(c_index,c_index+1) = tmp2;
   C33(c_index+1,c_index) = rxHalf;
   C33(c_index+1,c_index+1) = rxHalf;

   for ix = 1:length(is)
       m = 3*(is(ix) + Nx*it(ix));
       idex = idex+1;
       lin(idex) = n3;
       col(idex) = n3 + m;
       ele(idex) = C33(c_index+is(ix),c_index+it(ix));
   end

   
   for j = 3:Nx-1 % k = 1
       
       Kw = sqrt(-Y(k,j)*Z(k,j)); % on veut grandeur de K
       tmp1 = Kw*Kw - ky2Half + 2.0*1i*Kw/delta_z - rx;
       tmp2 = Kw*Kw - ky2Half - 2.0*1i*Kw/delta_z - rx;
       
       n1 = 3*(j + Nx*(k-1))-2;    % ligne assoiée à x
       n3 = 3*(j + Nx*(k-1));      % ligne assoicée à Jy
       is = [-1 -1 0 0 1 1];
       it = [0 1 0 1 0 1];
       
        %mêmes équations pour x et y
        
       C11(c_index-1,c_index) = rxHalf;
       C11(c_index-1,c_index+1)= rxHalf;
       C11(c_index,c_index)= tmp1;
       
       C11(c_index,c_index+1) = tmp2;
       C11(c_index+1,c_index) = rxHalf;
       C11(c_index+1,c_index+1) = rxHalf;
       
       vect_lin = [n1 n3];
       for ix = 1:length(is)
           for idex_line = 1:2
               n = vect_lin(idex_line);
               m = 3*(is(ix) + Nx*it(ix));
               idex = idex+1;
               lin(idex) = n;
               col(idex) = n + m;
               ele(idex) = C11(c_index+is(ix),c_index+it(ix));
           end
       end
   end
    
    
   %%%%%%%%%%%%%%%%% FIN FRONTIÈRE AU BAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%% FRONTIÈRE À GAUCHE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Pour composantes Jy et Jz, on passe par Engquist-Majda
    % Jx n'est pas couplé au système à cet endroit, donc aucune équation
    % associée...
    
    j = 1; 
    
    for k = 2:Nz-2
        Kw = sqrt(-Y(k,j)*Z(k,j)); % on veut grandeur de K
        tmp1 = Kw*Kw - ky2Half + 2.0*1i*Kw/delta_x - rz;
        tmp2 = Kw*Kw - ky2Half - 2.0*1i*Kw/delta_x - rz;
        
        n2 = 3*(j + Nx*(k-1))-1;    % ligne assoicée à z
        n3 = 3*(j + Nx*(k-1));      % ligne assoicée à y
        is = [0 0 0 1 1 1];
        it = [-1 0 1 -1 0 1];
        
        %mêmes équations pour z et y
        
        C22(c_index,c_index-1) = rzHalf;
        C22(c_index,c_index)= tmp1;
        C22(c_index,c_index+1)=rzHalf;
        
        C22(c_index+1,c_index-1) = rzHalf;
        C22(c_index+1,c_index) = tmp2;
        C22(c_index+1,c_index+1) = rzHalf;
       
        vect_lin = [n2 n3];
        for ix = 1:length(is)
            for idex_line = 1:2
                n = vect_lin(idex_line);
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n;
                col(idex) = n + m;
                ele(idex) = C22(c_index+is(ix),c_index+it(ix));
            end
        end
    end
    
    k = Nz-1; % j toujours = à 1
    
    Kw = sqrt(-Y(k,j)*Z(k,j)); % on veut grandeur de K
    tmp1 = Kw*Kw - ky2Half + 2.0*1i*Kw/delta_x - rz;
    tmp2 = Kw*Kw - ky2Half - 2.0*1i*Kw/delta_x - rz;
    
    n3 = 3*(j + Nx*(k-1));      % ligne assoicée à Jy; Jz traité comme un coin
    is = [0 0 0 1 1 1];
    it = [-1 0 1 -1 0 1];
    
    %mêmes équations pour z et y
    
    C33(c_index,c_index-1) = rzHalf;
    C33(c_index,c_index)= tmp1;
    C33(c_index,c_index+1)=rzHalf;
    
    C33(c_index+1,c_index-1) = rzHalf;
    C33(c_index+1,c_index) = tmp2;
    C33(c_index+1,c_index+1) = rzHalf;
  
    for ix = 1:length(is)
        m = 3*(is(ix) + Nx*it(ix));
        idex = idex+1;
        lin(idex) = n3;
        col(idex) = n3 + m;
        ele(idex) = C33(c_index+is(ix),c_index+it(ix));
    end
    
   %%%%%%%%%%%%%%% FIN FRONTIÈRE À GAUCHE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%% FRONTIÈRE À DROITE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    j = Nx;% frontière à droite
       
    % Pour composantes Jy et Jz, on passe par Engquist-Majda
    % Pour Jx, on traite comme un point intérieur
    for k = 2:Nz-2
        Kw = sqrt(-Y(k,j)*Z(k,j)); % on veut grandeur de K
        tmp1 = -Kw*Kw + ky2Half + 2.0*1i*Kw/delta_x + rz;
        tmp2 = -Kw*Kw + ky2Half - 2.0*1i*Kw/delta_x + rz;
        
        n2 = 3*(j + Nx*(k-1))-1;    % ligne assoicée à z
        n3 = 3*(j + Nx*(k-1));      % ligne assoicée à y
        is = [-1 -1 -1 0 0 0];
        it = [-1 0 1 -1 0 1];
        
        %mêmes équations pour z et y
        C22(c_index-1,c_index-1) = -rzHalf;
        C22(c_index-1,c_index) = tmp1;
        C22(c_index-1,c_index+1) = -rzHalf;
        
        C22(c_index,c_index-1) = -rzHalf;
        C22(c_index,c_index)= tmp2;
        C22(c_index,c_index+1)=-rzHalf;
       
        vect_lin = [n2 n3];
        for ix = 1:length(is)
            for idex_line = 1:2
                n = vect_lin(idex_line);
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n;
                col(idex) = n + m;
                ele(idex) = C22(c_index+is(ix),c_index+it(ix));
            end
        end
    end
 
    k = Nz-1; %j encore = à Nx...
    Kw = sqrt(-Y(k,j)*Z(k,j)); % on veut grandeur de K
    tmp1 = -Kw*Kw + ky2Half + 2.0*1i*Kw/delta_x + rz;
    tmp2 = -Kw*Kw + ky2Half - 2.0*1i*Kw/delta_x + rz;
  
    n3 = 3*(j + Nx*(k-1));      % ligne assoicée à y, on traite z omme un coin
    is = [-1 -1 -1 0 0 0];
    it = [-1 0 1 -1 0 1];
    
    C33(c_index-1,c_index-1) = -rzHalf;
    C33(c_index-1,c_index) = tmp1;
    C33(c_index-1,c_index+1) = -rzHalf;
    C33(c_index,c_index-1) = -rzHalf;
    C33(c_index,c_index)= tmp2;
    C33(c_index,c_index-1)=-rzHalf;
   
    for ix = 1:length(is)   
        m = 3*(is(ix) + Nx*it(ix));
        idex = idex+1;
        lin(idex) = n3;
        col(idex) = n3 + m;
        ele(idex) = C33(c_index+is(ix),c_index+it(ix));
    end
 
    %%%%%%%%%%%%%%% FIN FRONTIÈRE À DROITE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%%     Les points intérieurs
    %%% 
    %%% Pour ces points, on peut utiliser la discrétisation proposée
    %%% aux équations A6a, A6b, A6c de l'article d'Ellefsen. Ces équations 
    %%% sont codées par  9 fonction f_C11 à f_C33 un peu plus bas
    %%%
    %%% Cas particuliers:
    %%%
    %%%     Lorsque j = Nx et qu'il ne s'agit pas d'une arête, l'équation 
    %%%     associée à Jx peut se traiter comme un point intérieur.
    %%%
    %%%     Même chose pour k = 1 et l'équation en Jz. 
    %%%     
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    % on traite les deux cas particuliers en premier...
    
    j = Nx; % frontière à droite, équation pour Jx
    
    for k = 2:(Nz-1)
        % Pour composante X, on peut prendre les fonctions f_C11
        % f_C12, f_C13
        n1 = 3*(j + Nx*(k-1))-2; % la ligne à remplir pour X
        C11 = f_C11(Y,j,k,Ky,Zinv,XIz,delta_z);
        C12 = f_C12(j,k,Zinv,XIx,XIz,delta_x,delta_z);
        C13 = f_C13(j,k,Ky,Zinv,XIx,delta_x);
        
        is = [0 0 0];
        it = [-1 0 1];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n1;
            col(idex) = n1 + m;
            ele(idex) =  C11(c_index+is(ix),c_index+it(ix));
        end
        is = [-1 -1 0 0];
        it = [-1 0 -1 0];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n1;
            col(idex) = n1 + m + 1;
            ele(idex) = C12(c_index+is(ix),c_index+it(ix));
        end
        is = [-1 0];
        it = [0 0];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n1;
            col(idex) = n1 + m + 2;
            ele(idex) = C13(c_index+is(ix),c_index+it(ix));
        end
    end
    
    k = 1; % frontière en bas, équation pour Jz
    
    for j = 2:(Nx-1)
        n2 = 3*(j + Nx*(k-1))-1; % ligne pour Jz
        C21 = f_C21(j,k,Zinv,XIx,XIz,delta_x,delta_z);
        C22 = f_C22(Y,j,k,Ky,Zinv,XIx,delta_x);
        C23 = f_C23(j,k,Ky,Zinv,XIz,delta_z);
        
        is = [1 1 0 0];
        it = [1 0 1 0];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n2;
            col(idex) = n2 + m - 1;
            ele(idex) = C21(c_index+is(ix),c_index+it(ix));
        end
        is = [-1 0 1];
        it = [0 0 0];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n2;
            col(idex) = n2 + m;
            ele(idex) = C22(c_index+is(ix),c_index+it(ix));
        end
        is = [0 0];
        it = [1 0];
        
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n2;
            col(idex) = n2 + m + 1;
            ele(idex) = C23(c_index+is(ix),c_index+it(ix));
        end
        
    end
    
    
    %%% Patch à repenser...
    %
    % On ne peut pas utiliser les coefficients optimaux pour les points de 
    % l'avant dernière rangée ou l'avant dernière colonne
    %
    %
    jj = 2:Nx-1;
    kk1 = [2:Nx-1]./[2:Nx-1]*(Nz-1); %#ok<*NBRAK>
    kk2 = [2:Nx-1]./[2:Nx-1]*(Nz-2);
    
    temp = 2:Nz-3;
    %jj = [jj temp./temp*2 temp./temp*(Nx-1)];
    %kk = [kk temp temp];
    
    jj = [jj jj temp./temp*(Nx-2) temp./temp*(Nx-1)];
    kk = [kk1 kk2 temp temp];
   
    
    for ii = 1:length(jj)
        j = jj(ii);
        k = kk(ii);
        
        n1 = 3*(j + Nx*(k-1))-2; % les trois lignes à remplir pour j et k
        n2 = 3*(j + Nx*(k-1))-1;
        n3 = 3*(j + Nx*(k-1));
        
        C11 = f_C11(Y,j,k,Ky,Zinv,XIz,delta_z);
        C12 = f_C12(j,k,Zinv,XIx,XIz,delta_x,delta_z);
        C13 = f_C13(j,k,Ky,Zinv,XIx,delta_x);
        C21 = f_C21(j,k,Zinv,XIx,XIz,delta_x,delta_z);
        C22 = f_C22(Y,j,k,Ky,Zinv,XIx,delta_x);
        C23 = f_C23(j,k,Ky,Zinv,XIz,delta_z);
        C31 = f_C31(j,k,Ky,Zinv,XIx,delta_x);
        C32 = f_C32(j,k,Ky,Zinv,XIz,delta_z);
        C33 = f_C33(Y,j,k,Zinv,XIx,XIz,delta_x,delta_z);
        
        
        % Les coefficients pour c0 (c11, c22 et c33)
        is = [0 0 0];
        it = [-1 0 1];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n1;
            col(idex) = n1 + m;
            ele(idex) =  C11(c_index+is(ix),c_index+it(ix));
        end
        
        is = [-1 0 1];
        it = [0 0 0];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n2;
            col(idex) = n2 + m;
            ele(idex) = C22(c_index+is(ix),c_index+it(ix));
        end
        
        is = [0 -1 0 1 0];
        it = [0 0 -1 0 1];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n3;
            col(idex) = n3 + m;
            ele(idex) = C33(c_index+is(ix),c_index+it(ix));
        end
        
        % Les coefficients pour c1 (c12 et c21)
        is = [-1 -1 0 0];
        it = [-1 0 -1 0];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n1;
            col(idex) = n1 + m + 1;
            ele(idex) = C12(c_index+is(ix),c_index+it(ix));
        end
        
        is = [1 1 0 0];
        it = [1 0 1 0];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n2;
            col(idex) = n2 + m - 1;
            ele(idex) = C21(c_index+is(ix),c_index+it(ix));
        end
        
        % Les coefficients pour c13 et c31
        is = [-1 0];
        it = [0 0];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n1;
            col(idex) = n1 + m + 2;
            ele(idex) = C13(c_index+is(ix),c_index+it(ix));
        end
        
        is = [1 0];
        it = [0 0];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n3;
            col(idex) = n3 + m - 2;
            ele(idex) = C31(c_index+is(ix),c_index+it(ix));
        end
        
        % Les coefficients pour  c23 et c32
        is = [0 0];
        it = [1 0];
        
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n2;
            col(idex) = n2 + m + 1;
            ele(idex) = C23(c_index+is(ix),c_index+it(ix));
        end
        
        is = [0 0];
        it = [-1 0];
        for ix = 1:length(is)
            m = 3*(is(ix) + Nx*it(ix));
            idex = idex+1;
            lin(idex) = n3;
            col(idex) = n3 + m - 1;
            ele(idex) = C32(c_index+is(ix),c_index+it(ix));
        end
        
    end
    
    
    %%%% Tous les autres points intérieurs... 
    
    for k = 2:(Nz-3)
        for j = 2:(Nx-3)
            n1 = 3*(j + Nx*(k-1))-2; % les trois lignes à remplir pour j et k
            n2 = 3*(j + Nx*(k-1))-1;
            n3 = 3*(j + Nx*(k-1));
            
            C11 = f_C11Opt(Y,j,k,Ky,Zinv,XIz,delta_z,beta);
  %          C11 = f_C11(Y,j,k,Ky,Zinv,XIz,delta_z);
            C12 = f_C12(j,k,Zinv,XIx,XIz,delta_x,delta_z);
            C13 = f_C13(j,k,Ky,Zinv,XIx,delta_x);
            C21 = f_C21(j,k,Zinv,XIx,XIz,delta_x,delta_z);
            C22 = f_C22Opt(Y,j,k,Ky,Zinv,XIx,delta_x,beta);
            
       %     C22 = f_C22(Y,j,k,Ky,Zinv,XIx,delta_x);
            C23 = f_C23(j,k,Ky,Zinv,XIz,delta_z);
            C31 = f_C31(j,k,Ky,Zinv,XIx,delta_x);
            C32 = f_C32(j,k,Ky,Zinv,XIz,delta_z);
            C33 = f_C33Opt(Y,j,k,Zinv,XIx,XIz,delta_x,delta_z,beta);
           % C33 = f_C33(Y,j,k,Zinv,XIx,XIz,delta_x,delta_z);
            
            % Les coefficients pour c0 (c11, c22 et c33)
            is = [0 0 0 1 1 1 -1 -1 -1];
            it = [-1 0 1 -1 0 1 -1 0 1];
            for ix = 1:length(is)
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n1;
                col(idex) = n1 + m;
                ele(idex) =  C11(c_index+is(ix),c_index+it(ix));
            end
            
            is = [0 0 0 1 1 1 -1 -1 -1];
            it = [-1 0 1 -1 0 1 -1 0 1];
            for ix = 1:length(is)
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n2;
                col(idex) = n2 + m;
                ele(idex) = C22(c_index+is(ix),c_index+it(ix));
            end
            
            is = [0 0 0 1 1 1 -1 -1 -1];
            it = [-1 0 1 -1 0 1 -1 0 1];
            for ix = 1:length(is)
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n3;
                col(idex) = n3 + m;
                ele(idex) = C33(c_index+is(ix),c_index+it(ix));
            end
            
            % Les coefficients pour c1 (c12 et c21)
            is = [-1 -1 0 0];
            it = [-1 0 -1 0];
            for ix = 1:length(is)
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n1;
                col(idex) = n1 + m + 1;
                ele(idex) = C12(c_index+is(ix),c_index+it(ix));
            end
            
            is = [1 1 0 0];
            it = [1 0 1 0];
            for ix = 1:length(is)
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n2;
                col(idex) = n2 + m - 1;
                ele(idex) = C21(c_index+is(ix),c_index+it(ix));
            end
            
            % Les coefficients pour c13 et c31
            is = [-1 0];
            it = [0 0];
            for ix = 1:length(is)
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n1;
                col(idex) = n1 + m + 2;
                ele(idex) = C13(c_index+is(ix),c_index+it(ix));
            end
            
            is = [1 0];
            it = [0 0];
            for ix = 1:length(is)
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n3;
                col(idex) = n3 + m - 2;
                ele(idex) = C31(c_index+is(ix),c_index+it(ix));
            end
            
            % Les coefficients pour  c23 et c32
            is = [0 0];
            it = [1 0];
            
            for ix = 1:length(is)
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n2;
                col(idex) = n2 + m + 1;
                ele(idex) = C23(c_index+is(ix),c_index+it(ix));
            end
            
            is = [0 0];
            it = [-1 0];
            for ix = 1:length(is)
                m = 3*(is(ix) + Nx*it(ix));
                idex = idex+1;
                lin(idex) = n3;
                col(idex) = n3 + m - 1;
                ele(idex) = C32(c_index+is(ix),c_index+it(ix));
            end
        end
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Les fonctions pour points particuliers
 
 function K = f_C11(Y,j,k,Ky,Zinv,XIz,delta_z)
    % Coefficients sur la diagonale C11 et C22
 
    K = zeros(3,3)*(1 + 1i);

    c_index = 2; %
 
    Zinv_0_Pd = 1/2*(Zinv(k,j) + Zinv(k+1,j));
    Zinv_0_Md = 1/2*(Zinv(k-1,j) + Zinv(k,j));
  
    K(c_index,c_index) = Y(k,j) + (Zinv_0_Pd/XIz(k,j,3)+Zinv_0_Md/XIz(k,j,1))/XIz(k,j,2)/delta_z^2 + Ky^2*Zinv(k,j); %Ex_0,0
    K(c_index,c_index + 1) = -Zinv_0_Pd/delta_z^2/XIz(k,j,2)/XIz(k,j,3); %Ex_0,+1
    K(c_index,c_index - 1) = -Zinv_0_Md/delta_z^2/XIz(k,j,2)/XIz(k,j,1); %Ex_0,-1
 end
    
 function K = f_C12(j,k,Zinv,XIx,XIz,delta_x,delta_z)
    % Coefficients sur la diagonale C11 et C22
 
    K = zeros(3,3)*(1 + 1i);

    c_index = 2; %
    
    Zinv_0_Pd = 1/2*(Zinv(k,j) + Zinv(k+1,j));
    Zinv_0_Md = 1/2*(Zinv(k,j) + Zinv(k-1,j));
    XIx_0_Pd = 1/2*(XIx(k,j,2) + XIx(k+1,j,2));
    XIx_0_Md = 1/2*(XIx(k,j,2) + XIx(k-1,j,2));
    
    K(c_index,c_index) = Zinv_0_Pd/delta_z/delta_x/XIx_0_Pd/XIz(k,j,2); %E_z0x0
    K(c_index,c_index - 1) = -Zinv_0_Md/delta_z/delta_x/XIx_0_Md/XIz(k,j,2); %E_z0,-1
    K(c_index-1,c_index) = -Zinv_0_Pd/delta_z/delta_x/XIx_0_Pd/XIz(k,j,2); %E_z-1,0
    K(c_index-1,c_index-1) = Zinv_0_Md/delta_z/delta_x/XIx_0_Md/XIz(k,j,2); %E_z-1,-1
 end
 
 function K = f_C13(j,k,Ky,Zinv,XIx,delta_x)
    % Coefficients sur la diagonale C11 et C22
 
    K = zeros(3,3)*(1 + 1i);
    c_index = 2; %
   
    K(c_index,c_index) = 1i*Ky*Zinv(k,j)/delta_x/XIx(k,j,2); %E_y0,0
    K(c_index-1,c_index) = -1i*Ky*Zinv(k,j)/delta_x/XIx(k,j,2); %E_y-1,0
 end
 
 function K = f_C21(j,k,Zinv,XIx,XIz,delta_x,delta_z)
    % Coefficients sur la diagonale C11 et C22
 
    K = zeros(3,3)*(1 + 1i);
    c_index = 2; %
 
    XIx_d_d = 1/2*(XIx(k,j,3) + XIx(k+1,j,3));
    Zinv_1_d = 1/2*(Zinv(k,j+1) + Zinv(k+1,j+1));
    Zinv_0_d = 1/2*(Zinv(k,j) + Zinv(k+1,j));
    
    K(c_index,c_index) = Zinv_0_d/delta_x/delta_z/XIx_d_d/XIz(k,j,3); %E_x0,0
    K(c_index,c_index + 1) = -Zinv_0_d/delta_x/delta_z/XIx_d_d/XIz(k,j,3); %E_x0,+1
    K(c_index+1,c_index) = -Zinv_1_d/delta_x/delta_z/XIx_d_d/XIz(k,j+1,3); %E_x+1,0
    K(c_index+1,c_index+1) = Zinv_1_d/delta_x/delta_z/XIx_d_d/XIz(k,j+1,3); %E_x+1,+1
 end
    
 function K = f_C22(Y,j,k,Ky,Zinv,XIx,delta_x)
    % Coefficients sur la diagonale C11 et C22
 
    K = zeros(3,3)*(1 + 1i);

    c_index = 2; %

    D = 1/4*(Zinv(k+1,j+1) + Zinv(k,j+1) + Zinv(k+1,j) + Zinv(k,j));
    
    XIx_d_d = 1/2*(XIx(k,j,3) + XIx(k+1,j,3));
    XIx_1_d = 1/2*(XIx(k,j,4) + XIx(k+1,j,4));
    XIx_0_d = 1/2*(XIx(k,j,2) + XIx(k+1,j,2));
    Zinv_1_d = 1/2*(Zinv(k,j+1) + Zinv(k+1,j+1));
    Zinv_0_d = 1/2*(Zinv(k,j) + Zinv(k+1,j));
  
    K(c_index,c_index) = 1/4*(Y(k+1,j+1) + Y(k,j+1) + Y(k+1,j) + Y(k,j)) + Ky^2*D + (Zinv_0_d/XIx_0_d + Zinv_1_d/XIx_1_d)/XIx_d_d/delta_x^2; %E_z00
    K(c_index-1,c_index) = -Zinv_0_d/delta_x^2/XIx_d_d/XIx_0_d; %E_z-1,0
    K(c_index+1,c_index) = -Zinv_1_d/delta_x^2/XIx_d_d/XIx_1_d; %E_z+1,0
 end
 
 function K = f_C23(j,k,Ky,Zinv,XIz,delta_z)
    % Coefficients sur la diagonale C11 et C22
 
    K = zeros(3,3)*(1 + 1i);
    c_index = 2; %

    Dz = 1/4*(Zinv(k+1,j+1) + Zinv(k,j+1) + Zinv(k+1,j) + Zinv(k,j));
   
    K(c_index,c_index) = -1i*Ky*Dz/delta_z/XIz(k,j,3); %E_y0,0
    K(c_index,c_index+1) = 1i*Ky*Dz/delta_z/XIz(k,j,3); %E_y0,+1
 end
 
  function K = f_C31(j,k,Ky,Zinv,XIx,delta_x)
    % Coefficients sur la diagonale C11 et C22
 
    K = zeros(3,3)*(1 + 1i);

    c_index = 2; %
    
    K(c_index,c_index) = -1i*Ky*Zinv(k,j)/delta_x/XIx(k,j,3); %E_x0,0
    K(c_index+1,c_index) = 1i*Ky*Zinv(k,j+1)/delta_x/XIx(k,j,3); %E_x+1,0
 end
    
 function K = f_C32(j,k,Ky,Zinv,XIz,delta_z)
    % Coefficients sur la diagonale C11 et C22
 
    K = zeros(3,3)*(1 + 1i);
    c_index = 2; %
 
    D = 1/4*(Zinv(k+1,j+1) + Zinv(k,j+1) + Zinv(k+1,j) + Zinv(k,j));
    F = 1/4*(Zinv(k,j+1) + Zinv(k-1,j+1) + Zinv(k,j) + Zinv(k-1,j));
      
    K(c_index,c_index) = 1i*Ky*D/delta_z/XIz(k,j,2); %E_z00
    K(c_index,c_index-1) = -1i*Ky*F/delta_z/XIz(k,j,2); %E_z0,-1
 end
 
 function K = f_C33(Y,j,k,Zinv,XIx,XIz,delta_x,delta_z)
    % Coefficients sur la diagonale C11 et C22
 
    
    K = zeros(3,3)*(1 + 1i);
    c_index = 2; %
 
    XIz_d_0 = 1/2*(XIz(k,j,2) + XIz(k,j+1,2));
    XIz_d_Pd = 1/2*(XIz(k,j,3) + XIz(k,j+1,3));
    XIz_d_Md = 1/2*(XIz(k,j,1) + XIz(k,j+1,1));
    
    Zinv_d_Pd = 1/4*(Zinv(k+1,j+1) + Zinv(k,j+1) + Zinv(k+1,j) + Zinv(k,j));
    Zinv_d_Md = 1/4*(Zinv(k,j+1) + Zinv(k-1,j+1) + Zinv(k,j) + Zinv(k-1,j));
 
    K(c_index,c_index) = 1/2*(Y(k,j+1)+Y(k,j))+(Zinv_d_Pd/XIz_d_Pd + Zinv_d_Md/XIz_d_Md)/delta_z^2/XIz_d_0 + (Zinv(k,j)/XIx(k,j,2)+Zinv(k,j+1)/XIx(k,j,4))/delta_x^2/XIx(k,j,3); %E_y0,0
    K(c_index-1,c_index) = -Zinv(k,j)/delta_x^2/XIx(k,j,2)/XIx(k,j,3); %E_y-1,0
    K(c_index,c_index-1) = -Zinv_d_Md/delta_z^2/XIz_d_0/XIz_d_Md; %E_y0,-1
    K(c_index,c_index+1) = -Zinv_d_Pd/delta_z^2/XIz_d_0/XIz_d_Pd; %E_y0,+1
    K(c_index+1,c_index) = -Zinv(k,j+1)/delta_x^2/XIx(k,j,4)/XIx(k,j,3); %E_y+1,0
                            
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Les fonctions pour coeffients optimaux


 function K = f_C11Opt(Y,j,k,Ky,Zinv,XIz,delta_z,beta)
   
    K = zeros(3,3)*(1 + 1i);

    c_index = 2; %
    
    b = beta(1);
    c = beta(2);
    d = (1-beta(2))/2;
 
    Az0 = c*0.5*(Zinv(k,j) + Zinv(k+1,j))/XIz(k,j,3)/XIz(k,j,2)/delta_z^2;
    AzP1 = d*0.5*(Zinv(k,j+1) + Zinv(k+1,j+1))/XIz(k,j+1,3)/XIz(k,j,2)/delta_z^2;
    AzM1 = d*0.5*(Zinv(k,j-1) + Zinv(k+1,j-1))/XIz(k,j-1,3)/XIz(k,j,2)/delta_z^2;
    
    Bz0 = c*0.5*(Zinv(k-1,j) + Zinv(k,j))/XIz(k,j,1)/XIz(k,j,2)/delta_z^2;
    BzP1 = d*0.5*(Zinv(k-1,j+1) + Zinv(k,j+1))/XIz(k,j+1,1)/XIz(k,j,2)/delta_z^2;
    BzM1 = d*0.5*(Zinv(k-1,j-1) + Zinv(k,j-1))/XIz(k,j-1,1)/XIz(k,j,2)/delta_z^2;
    
    K(c_index,c_index) = Ky^2*Zinv(k,j) + b*Y(k,j) + Az0 + Bz0; %Ex_0,0 
    K(c_index+1,c_index) = (1-b)/4*Y(k,j+1) + AzP1 + BzP1;  %Ex_1,0
    K(c_index-1,c_index) = (1-b)/4*Y(k,j-1) + AzM1 + BzM1;  %Ex_-1,0
    
                         
    K(c_index,c_index + 1) = (1-b)/4*Y(k+1,j)- Az0; %Ex_0,+1
    K(c_index+1,c_index + 1) = -AzP1; %Ex_+1,+1
    K(c_index-1,c_index + 1) = -AzM1; %Ex_-1,+1

            
    K(c_index,c_index - 1) = (1-b)/4*Y(k-1,j)-Bz0; %Ex_0,-1 
    K(c_index+1,c_index - 1) = -BzP1; %Ex_+1,-1 
    K(c_index-1,c_index - 1) = -BzM1; %Ex_-1,-1 
                                         
 end
    
 
    
 function K = f_C22Opt(Y,j,k,Ky,Zinv,XIx,delta_x,beta)
    % Coefficients sur la diagonale C11 et C22
 
    K = zeros(3,3)*(1 + 1i);

    c_index = 2; %
    
    b = beta(1);
    c = beta(2);
    d = (1-beta(2))/2;
    
    D = 0.25*(Zinv(k+1,j+1) + Zinv(k,j+1) + Zinv(k+1,j) + Zinv(k,j));
    
    XIx_d_d = 1/2*(XIx(k,j,3) + XIx(k+1,j,3));
    XIx_1_d = 1/2*(XIx(k,j,4) + XIx(k+1,j,4));
    XIx_1_dP1 = 1/2*(XIx(k+1,j,4) + XIx(k+2,j,4));
    XIx_1_dM1 = 1/2*(XIx(k,j,4) + XIx(k-1,j,4));
    
    XIx_0_d = 1/2*(XIx(k,j,2) + XIx(k+1,j,2));
    XIx_0_dP1 = 1/2*(XIx(k+1,j,2) + XIx(k+2,j,2));
    XIx_0_dM1 = 1/2*(XIx(k,j,2) + XIx(k-1,j,2));
    
    Zinv_1_d = 1/2*(Zinv(k,j+1) + Zinv(k+1,j+1));
    Zinv_1_dP1 = 1/2*(Zinv(k+1,j+1) + Zinv(k+2,j+1));
    Zinv_1_dM1 = 1/2*(Zinv(k,j+1) + Zinv(k-1,j+1));
    Zinv_0_d = 1/2*(Zinv(k,j) + Zinv(k+1,j));
    Zinv_0_dP1 = 1/2*(Zinv(k+1,j) + Zinv(k+2,j));
    Zinv_0_dM1 = 1/2*(Zinv(k,j) + Zinv(k-1,j));
    
    Ax0 = c*Zinv_1_d/XIx_1_d/XIx_d_d/delta_x^2;
    AxP1 = d*Zinv_1_dP1/XIx_1_dP1/XIx_d_d/delta_x^2;
    AxM1 = d*Zinv_1_dM1/XIx_1_dM1/XIx_d_d/delta_x^2;
    
    Bx0 = c*Zinv_0_d/XIx_0_d/XIx_d_d/delta_x^2;
    BxP1 = d*Zinv_0_dP1/XIx_0_dP1/XIx_d_d/delta_x^2;
    BxM1 = d*Zinv_0_dM1/XIx_0_dM1/XIx_d_d/delta_x^2;
     
    Y_h_h = 0.25*(Y(k+1,j+1) + Y(k,j+1) + Y(k+1,j+1) + Y(k,j));
    Y_hP1_h = 0.25*(Y(k+1,j+2) + Y(k,j+2) + Y(k+1,j+1) + Y(k,j+1));
    Y_hM1_h = 0.25*(Y(k+1,j) + Y(k,j) + Y(k+1,j-1) + Y(k,j-1));
    Y_h_hP1 = 0.25*(Y(k+2,j+1) + Y(k+1,j+1) + Y(k+2,j) + Y(k+1,j));
    Y_h_hM1 = 0.25*(Y(k,j+1) + Y(k-1,j+1) + Y(k,j) + Y(k-1,j));
    
  
    K(c_index,c_index) = Ky^2*D + b*Y_h_h + Ax0 + Bx0;     %E_z00                  
    K(c_index,c_index+1) = (1-b)/4*Y_h_hP1 + AxP1 + BxP1;   %E_z0+1             
    K(c_index,c_index-1) = (1-b)/4*Y_h_hM1 + AxM1 + BxM1;   %E_z0-1   
                           
    K(c_index-1,c_index) = (1-b)/4*Y_hM1_h - Bx0;   %E_z-1,0
    K(c_index-1,c_index+1) = -BxP1;                %E_z-1,+1
    K(c_index-1,c_index-1) = -BxM1;                %E_z-1,-1
                                                                                               
    K(c_index+1,c_index) = (1-b)/4*Y_hP1_h - Ax0;   %E_z+1,0
    K(c_index+1,c_index+1) = -AxP1;                %E_z+1,+1
    K(c_index+1,c_index-1) = -AxM1;                %E_z+1,-1                    
      
 end
 

   
 function K = f_C33Opt(Y,j,k,Zinv,XIx,XIz,delta_x,delta_z,beta)
    % Coefficients sur la diagonale C11 et C22
 
    b = beta(1);
    c = beta(2);
    d = (1-beta(2))/2;
    
    K = zeros(3,3)*(1 + 1i);
    c_index = 2; %
 
    XIz_d_0 = 1/2*(XIz(k,j,2) + XIz(k,j+1,2));
    XIz_d_Pd = 1/2*(XIz(k,j,3) + XIz(k,j+1,3));
    XIz_dP1_Pd = 1/2*(XIz(k,j+1,3) + XIz(k,j+2,3));
    XIz_dM1_Pd = 1/2*(XIz(k,j,3) + XIz(k,j-1,3));
    
    XIz_d_Md = 1/2*(XIz(k,j,1) + XIz(k,j+1,1));
    XIz_dP1_Md = 1/2*(XIz(k,j+1,1) + XIz(k,j+2,1)); 
    XIz_dM1_Md = 1/2*(XIz(k,j,1) + XIz(k,j-1,1));
    
    Zinv_d_Pd = 1/4*(Zinv(k+1,j+1) + Zinv(k,j+1) + Zinv(k+1,j) + Zinv(k,j));
    Zinv_dP1_Pd = 1/4*(Zinv(k+1,j+1) + Zinv(k,j+1) + Zinv(k+1,j+2) + Zinv(k,j+2));
    Zinv_dM1_Pd = 1/4*(Zinv(k+1,j) + Zinv(k,j) + Zinv(k+1,j-1) + Zinv(k,j-1));
    
    Zinv_d_Md = 1/4*(Zinv(k,j+1) + Zinv(k-1,j+1) + Zinv(k,j) + Zinv(k-1,j));
    Zinv_dP1_Md = 1/4*(Zinv(k,j+1) + Zinv(k-1,j+1) + Zinv(k,j+2) + Zinv(k-1,j+2));
    Zinv_dM1_Md = 1/4*(Zinv(k,j) + Zinv(k-1,j) + Zinv(k,j-1) + Zinv(k-1,j-1));
    
    
    Ax0 = c*Zinv(k,j+1)/delta_x^2/XIx(k,j,4)/XIx(k,j,3);
    AxP1 = d*Zinv(k+1,j+1)/delta_x^2/XIx(k+1,j,4)/XIx(k,j,3);
    AxM1 = d*Zinv(k-1,j+1)/delta_x^2/XIx(k-1,j,4)/XIx(k,j,3);
    
    Bx0 = c*Zinv(k,j)/delta_x^2/XIx(k,j,2)/XIx(k,j,3);
    BxP1 = d*Zinv(k+1,j)/delta_x^2/XIx(k+1,j,2)/XIx(k,j,3);
    BxM1 = d*Zinv(k-1,j)/delta_x^2/XIx(k-1,j,2)/XIx(k,j,3);
    
    Az0 = c*Zinv_d_Pd/XIz_d_0/XIz_d_Pd/delta_z^2;
    AzP1 = d*Zinv_dP1_Pd/XIz_d_0/XIz_dP1_Pd/delta_z^2;
    AzM1 = d*Zinv_dM1_Pd/XIz_d_0/XIz_dM1_Pd/delta_z^2;
    
    Bz0 = c*Zinv_d_Md/XIz_d_0/XIz_d_Md/delta_z^2;
    BzP1 = d*Zinv_dP1_Md/XIz_d_0/XIz_dP1_Md/delta_z^2;
    BzM1 = d*Zinv_dM1_Md/XIz_d_0/XIz_dM1_Md/delta_z^2;
    
    
    Y_h_0 = 1/2*(Y(k,j+1)+Y(k,j));  
    Y_hP1_0 = 1/2*(Y(k,j+2)+Y(k,j+1));
    Y_hM1_0 = 1/2*(Y(k,j)+Y(k,j-1));
    Y_h_P1 = 1/2*(Y(k+1,j+1)+Y(k+1,j));
    Y_h_M1 = 1/2*(Y(k-1,j+1)+Y(k-1,j));
   
    K(c_index,c_index) = b*Y_h_0 + Ax0 + Bx0 + Az0 + Bz0; %E_y0,0                        
    K(c_index-1,c_index) = (1-b)/4*Y_hM1_0  + AzM1 + BzM1 - Bx0; %E_y-1,0                               
    K(c_index+1,c_index) =  (1-b)/4*Y_hP1_0 + AzP1 + BzP1 - Ax0; %E_y+1,0                             
    K(c_index,c_index-1) = (1-b)/4*Y_h_M1 - Bz0 + AxM1 + BxM1; %E_y0,-1
    K(c_index,c_index+1) = (1-b)/4*Y_h_P1 - Az0 + AxP1 + BxP1;  %E_y0,+1

    K(c_index+1,c_index+1) = -AxP1 - AzP1;
    
    K(c_index-1,c_index+1) = -BxP1 - AzM1;
    
    K(c_index+1,c_index-1) = -AxM1 - BzP1;
    
    K(c_index-1,c_index-1) = -BxM1 - BzM1;
  
 end
 