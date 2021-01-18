function [ M ] = trapz_matrix1(  X, varargin )
%% TRAPZ_MATRIX1, RECOVERING THE VECTOR FROM NUMERICAL INTEGRATION 
%     
%     M = trapz_matrix1(X) returns the underlying vector M
%     representing the numerical integration operator trapz over a vector X.
%
%     The output matrix M allows for converting a vectorized function
%     F(X) into the scalar G:
%           G(1,1) = M(1,1:Nx) * F(1:Nx,1)
%
%% FUNCTION OPTIONS
%
%     M = trapz_matrix1(X,V) returns the same matrix M, but with some
%     differential volume element V that depends on X. V is either a scalar
%     multiplied into each grid point equally, or is a vector that matches
%     the length of X.
%       
%     M = trapz_matrix1(X,METHOD) returns the same matrix M, but with
%     the method specified by METHOD. The default is trapz, but 'Simpsons'
%     is also supported for composite Simpson integration. (If Nx is even,
%     Simpsons 3/8 rule is used on the last three elements)
%     
%     Note: the order of METHOD and V does not matter. 
%     
    
%% Process inputs
    % Determine if V(vector) or METHOD (string) are passed
    
    % set default inputs 
    method = 1; % trapz
    V = []; % Empty 
    
    % Specify input type, if provided
    iarg = nargin;
    while (iarg > 1)
        
       % Test type
       if ( isstring(varargin{iarg-1}) || ...
            ischar(varargin{iarg-1}))
           
           % Extract Method
           method = varargin{iarg-1};
           
           % Test type
           if ( strcmpi( method, 'trapz' ) ) % skip
               method = 1;
           elseif (strcmpi( method, 'simpsons' ) )
               method = 2;
           else
               error('Bad string input for METHOD')
           end
           
       % Test vector
       elseif ( isvector(varargin{iarg-1}) )
           
           % Store Vmult Array
           V = varargin{iarg-1};
           
       % Error Case
       else
           error('Unknown VARARGIN')
       end
       
       % decrement iarg
       iarg = iarg - 1;
       
    end
    
%% Pre-process the Reference Data
    
    % Check if successful
    if (~isvector(X))
        error('X is not a vector')
    end
    
    % Reshape to vector, stored like target arrays (outer, in columns)
    xref = reshape( X, [], 1 );
    
    % Get lengths of the reference arrays, these can be different
    Nx = numel( xref );
    
    % Get cell width, note: cubic arrays must be uniform 
    dx(:,1) = diff( xref ); 
    
    % Tack on the last cell
    % ( this isn't used b/c points using it should be culled)
    dx( end + 1 ) = dx( end );
    
    % Make sure points don't overlap
    if any( dx==0 )
        error('grid not distinct')
    end
    
    % Make sure grids are sorted
    if ~issorted( xref )   
        error('grid is not sorted')
    end
    

   %%  Create Integration Coefficients V 
   
   % Initialize V if not provided
   if ( isempty(V) )
      V = ones(Nx,1);
   elseif ( isscalar(V) )
      V(1:Nx,1) = V(1);
   elseif ( ~isnumeric(V) )
       error('Provided V is not numeric')
   elseif ( length(V)~=Nx )
       error('Length of V does not match X')
   end 
   
   % Reshape
   V = reshape(V,[Nx,1]);
   
   % Multiply by axis delta
   V = V .* dx;
   
    % Fork Integration Method
    if (method == 1 ) % trapz

        % TrapZ Coefficients
        V(1  ) =   V(1  ) .*0.5;
        V(end) =   V(end) .*0.5;

    elseif (method == 2 ) % Simpsons

        % Set default end index for Simpsons
        Nv = length(V);
        
        % If even grid points (odd panels)
        if ( mod(Nv,2)==0 )
             
            % Reset limits for Simpson
             Nv  = Nv - 3; 
             
             % Save state of the last 4 terms
             Vend(1:4) = V(Nv:end);
             
        end
        
        % Simpsons Common Factor
        V(1:Nv) = V(1:Nv)./3.0;
        
        % Simpsons Multipliers
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
    
   %% Form of Final Integration Matrix
   M(1,1:Nx) = V;
   
end

