function [ M, Vout ] = trapz_matrix2(  X, Y, DIM, varargin )
%% TRAPZ_MATRIX2, RECOVERING THE MATRIX FROM NUMERICAL INTEGRATION 
%     
%     M = trapz_matrix2(X,Y,DIM) returns the underlying sparse matrix M
%     representing the numerical integration operator trapz over dimension
%     DIM on some 2D arbitrary coordinate map X/Y, where X and Y are either
%     vectors or matrices produced by ndgrid. Meshgrid is also accepted,
%     but is converted to ndgrid. DIM is a scalar of either 1 or 2,
%     signifying integration over either X or Y respectively.
%
%     **WARNING** Y IS NOT the Y typically seen in Trapz, (the solution
%     reference vector), but rather the geometric axis paired with X.
%
%     The output matrix M allows for converting a vectorized function
%     F(:)=vec(F(X,Y)) of size 1xNx*Ny, into G(X) or G(Y), depending on DIM
%     by:
%           G(1:Nx/Ny,1) = M(1:Nx/Ny,1:Ny*Nx) * F(1:Nx*Ny,1)
%     where DIM is the dimension which is collapsed, such that G=G(1:Ny,1)
%     is returned for DIM=1, and G=G(1:Nx,1) is returned for DIM=2.
%
%% FUNCTION OPTIONS
%
%     M = trapz_matrix2(X,Y,DIM,V) returns the same matrix M, but with some
%     differential volume element V that depends on the coordinate system
%     that X and Y represent. V is either a scalar applied to each term in
%     DIM, or is a vector that matches the length of the dimension DIM,
%     with implicit expansion across the other dimension. V does not include
%     the X-Axis spacing term (i.e. dx or dy)
%       
%     M = trapz_matrix2(X,Y,DIM,METHOD) returns the same matrix M, but with
%     the method specified by METHOD. The default is trapz, but 'Simpsons'
%     is also allowed for composite Simpson integration.
%
%     [ M, VOUT ] = trapz_matrix2(X,Y,DIM,METHOD) returns M and VOUT, where
%     M is split into M and VOUT, where M contains the same footprint but
%     only ones (a transformation matrix), and VOUT is the vector of
%     V=V(1:Nx*Ny,1) coefficients that can be multiplied elementwise with
%     F(1:Nx*Ny,1) prior to matrix mult with M. This is enabled to allow
%     access to the dimensional transformation matrix, which can be
%     inverted to map the reduced order G(1:Nz/Ny) back onto the full X/Y
%     space, if that is desired. The order of METHOD and V does not matter. 
%
%
%% Further Notes
%
%     **WARNING** trapz_matrix2 is for integrating along a SINGLE DIM of
%     X/Y. For 2D integration, trapz_matrix1 must be called, and the resulting
%     matrix N can be multiplied directly with matrix M to recover the full
%     integration or double contraction.
%     
%     **WARNING** Integration over Y is supported, but likely very
%     inefficient due to the stride. Best practice is to reshape to
%     integrate over X only. 
    
%% Process inputs
    % Determine if V(vector) or METHOD (string) are passed
    
    % set default inputs 
    method = 1; % trapz
    V = []; % Empty 
    
    % Specify input type, if provided
    iarg = nargin;
    while (iarg > 3)
        
       % Test type
       if ( isstring(varargin{iarg-3}) || ...
            ischar(varargin{iarg-3}))
           
           % Extract Method
           method = varargin{iarg-3};
           
           % Test type
           if ( strcmpi( method, 'trapz' ) ) % skip
               method = 1;
           elseif (strcmpi( method, 'simpsons' ) )
               method = 2;
           else
               error('Bad string input for METHOD')
           end
           
       % Test vector
       elseif ( isvector(varargin{iarg-3}) )
           
           % Store Vmult Array
           V = varargin{iarg-3};
           
       % Error Case
       else
           error('Unknown VARARGIN')
       end
       
       % decrement iarg
       iarg = iarg - 1;
       
    end
    
%% Pre-process the Reference Data
    
    % Get the correct 1D vectors from the given inputs
    [mtype, xref, yref ]= mesh_type( X, Y );
    
    % Check if successful
    if (~isvector(xref))
        error('X and Y are not of correct format')
    elseif ( mtype==3 )
        warning('Meshgrid format is not advised')
    end
    
    % Reshape to vector, stored like target arrays (outer, in columns)
    xref = reshape( xref, 1, [] );
    yref = reshape( yref, 1, [] );
    
    % Get lengths of the reference arrays, these can be different
    Nx = numel( xref );
    Ny = numel( yref );
    
    % Total length of a vectorized grid matrix
    %  (this is also the number of columns in M and rows in F)
    Nr = Nx*Ny;    
    
    % Get cell width, note: cubic arrays must be uniform 
    if (DIM==1) % integrate over X
        delta(:,1) = diff( xref ); 
        Nq = Ny; % Number of Rows in M matrix
        Nc = Nx; % Columns in M
    elseif (DIM==2) % over y
        delta(:,1) = diff( yref );    
        Nq = Nx; % Number of Rows in M matrix
        Nc = Ny; % Columns in M
    end
    
    % Tack on the last cell
    % ( this isn't used b/c points using it should be culled)
    delta( end + 1 ) = delta( end );
    
    % Test Scalar DIM
    if ( DIM~=1 && DIM ~=2)
        error('Bad DIM provided, must be one or two')
    end
    
    % Make sure points don't overlap
    if any( delta==0 )
        error('grid not distinct')
    end
    
    % Make sure grids are sorted
    if ~issorted( xref ) || ~issorted( yref )     
        error('grid is not sorted')
    end
    

   %%  Create Integration Coefficients V 
   
   % Initialize V if not provided
   if ( isempty(V) )
        V = ones(Nc,1);
   elseif ( isscalar(V) )
        V(1:Nc,1) = V(1,1);
   elseif ( ~isnumeric(V) )
       error('Provided V is not numeric')
   elseif ( length(V)~=Nc )
       error('Length of V does not match integration dimension')
   end 
   
   % Multiply by axis delta
   V = V .* delta;
   
    % Fork Method
    if (method == 1 ) % trapz

        % TrapZ Coefficients
        V(1  ) =   V(1  ) .*0.5;
        V(end) =   V(end) .*0.5;

    elseif (method == 2 ) % Simpsons

        % Set default end index for Simpsons
        Nv = length(V);
        
        % If even grid points (odd panels)
        if ( mod(Nv,2)==0 )
             Nv  = Nv - 3;
             Vend(1:4) = V(Nv:end);
        end
        
        % Common Factor
        V(1:Nv) = V(1:Nv)./3.0;
        
        % Multipliers
        V( 2:2:Nv-1  ) = 4.0 .* V( 2:2:Nv-1 ) ;
        V( 3:2:Nv-2  ) = 2.0 .* V( 3:2:Nv-2 ) ;
        
        % Finish Even case with Simpsons 3/8
        if ( Nv < length(V) )
            
            % Overlap term used in both stencils
            V(Nv)   = V(Nv)   + 0.375*Vend(1);
            
            % Remaining terms 
            V(Nv+1) = 1.1250*Vend(2);
            V(Nv+2) = 1.1250*Vend(3);
            V(Nv+3) =  0.375*Vend(4);
            
        end
        
    end
    
   %%  Create Integration Coefficients V 
   
    % Create Empty vector Vout
    Vout = ones(Nr,1);
    
    % Fork Axis 
    if (DIM==1) % integrate over X

        % Copy Arrays
        Vout = repmat(V,Ny,1);

        % Row Pointer
        I = reshape( repmat(1:Ny, Nx,1) , [] , 1 );
        
    % integrate over y
    elseif (DIM==2) 
            
        % Copy Arrays
        Vout = reshape( repmat(V,1,Nx)' , [],1 );
        
        % Row Pointer
        I = repmat([1:Nx]', Ny,1);
        
    end
    
%%  Create Transformation Matrix
    
    % Create the sparse array
    M = sparse( I, [1:Nr]', ones(Nr,1), Nq, Nr );
    
    % See if these should be merged or not
    if (nargout==1)
        M = bsxfun(@times, Vout' , M ); 
    end
    
end

