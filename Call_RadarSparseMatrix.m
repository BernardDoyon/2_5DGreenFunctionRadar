% ======================================================================
% Copyright (c) June 2019, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================
function Call_RadarSparseMatrix

    % This function is a dummy function to validate the building of the 
    % MEX file associated with RadarSparseMatrix.m 
    %
    %  To build the MEX file, you can use this dummy function to  
    %  automatically define input types with the MATLAB coder
    %  BUT you need to set the size of Y, Z, XIx, and XIz to 
    %  "undounded" (:inf) 
    %       Y : complex(double :inf x :inf)
    %       Z : complex(double :inf x :inf)
    %       XIx: complex(double :inf x :inf x 4)
    %       XIz: complex(double :inf x :inf x 4)
    
    Nx = 10;      
    Nz = 5;
    delta_x = 0.01;
    delta_z = 0.01;
    Ky = 1.0 + 1i/10;
    
    Y = ones(Nz,Nx)*(1 + 1i);
    Z = ones(Nz,Nx)*(1 + 1i);
    
    XIx = ones(Nz,Nx,4)*(1 + 1i);
    XIz = ones(Nz,Nx,4)*(1 + 1i); 
   
    [lin,col,ele] = RadarSparseMatrix(int32(Nx),int32(Nz),Ky,delta_x,...
                        delta_z,Y,Z,XIx,XIz);
    
end