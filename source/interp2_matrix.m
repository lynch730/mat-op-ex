function [ M, xout, yout, ob_pntr ] = interp2_matrix( X, Y, Xq, Yq, varargin )
%% INTERP2_MATRIX, RECOVERING THE MATRIX FROM INTERP2 
%     
%     M = interp2_matrix(X,Y,Xq,Yq) returns the underlying matrix M
%     representing the interp2 operator for a known X/Y grid and set of a
%     set of query points Xq/Yq. Matrix M is of size M(Nq,Nr), with Nq rows
%     of vectorized queries Xq(:)/Yq(:), and Nr columns of vectorized
%     reference grid points X(:)/Y(:). Out-of-Bounds queries are left
%     sparse in M (See ob_pntr).
%     
%     M can subsequently be used to obtain a vector P(Nq,1) of
%     Nq interpolation queries for a reference function F over X and Y:
%        P( 1:Nq,1 ) = M( 1:Nq, 1:Nr ) * F( 1:Nr,1 )
%     
%     While M could in theory be constructed to preserve the shapes of X/Y
%     and Xq/Yq in a high-order tensor that is doubly contracted over
%     F(X,Y), interp2_matrix vectorizes to maintain compatibility with
%     MATLAB matrix multiplication. It is left to the user to reshape M if
%     that is desired.
%     
%     WARNING: Because vectorization is in the default row-major order, X/Y
%     should either be vectors of equal length and uniform spacing, or
%     matrices generated from ndgrid. Meshgrid format is allowed for X/Y,
%     but is immediately converted with a warning. Likewise, any function
%     F(X,Y) created from meshgrid X/Y must be transposed BEFORE
%     vectorization in the above equation for P( 1:Nq,1 ). F(X,Y) generated
%     from ndgrid can simply be vectorized by F(:).
% 
%
%% FUNCTION OPTIONS
% 
%     M = interp2_matrix(X,Y,Xq,Yq,METHOD) specifies interpolation method,
%     which is linear by default. Other options are nearest and cubic.
%     Results should be identical to interp2, with the exception of cubic
%     cases where the stencil requires cells outside the grid. Matching
%     interp2 requires full calculation of not-a-knot spline, here we
%     approximate by extrapolating the missing cells with 2nd order
%     differences. Typical deviations from interp2 cubic are usually <1%.
%     
%       'nearest' - nearest neighbor interpolation 
%       'linear'  - bi-linear interpolation 
%       'cubic'   - bi-cubic convolution interpolation
%
%    [ M, XOUT, YOUT ] = interp2_matrix(...) returns the 1D X/Y position
%    arrays of each column in M. xout and yout are dimension Nr, are
%    equivalent to vectorizing an X/Y mesh created from ndgrid.
%    
%    [ M, XOUT, YOUT, OB_PNTR ] = interp2_matrix(...) returns an array of
%    the M matrix's row indicies corresponding to query points that are out
%    of bounds. These rows are left sparse. 
%
%% SOLUTION TECHNIQUE & SOURCE
%    
%    This solution is based heavily on the suggested technique in
%     <https://www.mathworks.com/matlabcentral/answers/573703-get-interpolation-transfer-relation-matrix-instead-of-interpolated-values>
%    which is expanded upon to include the cubic case 
%    
%    xref/yref are the 1D representation of X and Y. In terms of notation,
%    the x and y axes are indexed with i/j in the usual way. The stencil
%    "sten" are the grid points associated with individual interpolations.
%    Global "glbl" refers to indices/coordinates in the frame of the larger
%    xref/yref reference grid "ref".
%    
%    The general approach for each query is:
%       (1) Build the arrays i_glbl_pntr/j_glbl_pntr (integer arrays which
%           point each query to the nearest left-hand-side grid point in
%           xref/yref)
%       (2) Find the appropriate stencil (2x2 for linear and nearest,
%           4x4 cubic) of weights W, based on cell-relative position x/y
%       (3) Compute i_glbl_sten/j_glbl_sten which map each W to the
%           reference grid.
%       (4) Convert i_glbl_sten/j_glbl_sten into J, the vectorized position
%           in X/Y. Similar vectorize W from a 2D stencil.
%       (5) Construct the sparse matrix from W,I,J, where I is the row
%           index of M, the same for all stencil weights in each query.
%    The above (1-5) are all simultaneously performed for each query pair in
%    Xq/Yq, indexed on the outside of all of the above operations.
%    
%    Though not wildly slow, there is likely room to optimize the current
%    technique. However user accessibility is perhaps more important here.
%    The primary utility in obtaining the matrix is for minimizing the
%    overhead associated with repeated interp2 calls on a fixed grid, where
%    the M matrix need only be computed once. Note, this can be done to
%    some extent for conventional CPU operations by creating
%    P=griddedInterpolant (the underlying code interp2 wraps around), and
%    updating only the F.Values to avoid recomputing the weights.
%    Critically, this approach currently cannot be replicated for
%    gpuArrays, where such optimizations would be most useful. Of course,
%    interp2_matrix is still useful for analysis of Jacobians and other
%    applications that griddedInterpolant cannot reproduce.
%    
%    Note, some operations (like vectorized dx/dy and histcounts2), are
%    setup to enable future compatibility with unevenly spaced Cartesian
%    grids. This is trivial for bi-linear interpolation, but tricky and
%    potentially useless for bicubic

    %% Process Optional Arguments
    
    % Specify input type, if provided
    method = 'linear';
    if ( nargin > 4 )
       method = varargin{1};
    end
    
    % Read in custom BC's for cubic case
    cubic_bc_x = []; % default to extrapolation
    cubic_bc_y = [];
    if ( nargin == 7 )
       cubic_bc_x = varargin{2};
       cubic_bc_y = varargin{3};
    elseif (nargin == 6 ) % not enough supplied
        error('If applying cubic boundary conditions, need two arrays');
    end    
    
    %% Check the Interpolation Method & BC's
    
    % Fork interp string type
    if     strcmpi( method, 'nearest' )
       n_order = 0;
    elseif     strcmpi( method, 'linear'  )
       n_order = 1;
    elseif strcmpi( method, 'cubic'   )
       n_order = 2;
    else
        error('Bad Interp String, must be cubic or linear')
    end
    
    % Check to make sure cubic_bc's were correctly supplied
    if ( ~isempty(cubic_bc_x) && n_order~=2  )
        error('Cubic boundary conditions given when not using METHOD cubic')
        
    % Remaining tests if some BC's supplied
    elseif ( ~isempty(cubic_bc_x) && ~isempty(cubic_bc_y) )
        
        % Now test if both are not empty
        if ( isempty(cubic_bc_x) || isempty(cubic_bc_y)  )
            error('both boundary conditions are not supplied')
            
        % Test if dimensions and types are correct
        elseif( ~ismatrix(cubic_bc_x) )
            error('Cubic bc not a matrix of integers of order=2') 
            
        % Check outter sizes
        elseif ( size(cubic_bc_x,2) ~= 4  || ...
                 size(cubic_bc_y,2) ~= 4 )
            error('I/J not of size Nx4')
            
        end
        
    else % ensure no BC's are used
        cubic_bc_x = [];
        cubic_bc_y = [];
    end
    
    %% Pre-process the Query Data
    
    % Reshape to linear, with queries stored along outer index (columns)
    Xq = reshape( Xq, 1, [] );
    Yq = reshape( Yq, 1, [] );
    
    % Get the length of the linear arrays
    Nq = numel( Yq );
    
    % Test if lengths are the same
    if ( numel( Xq ) ~= Nq )
        error('Target X/Y arrays are not equal length')
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
    nx_ref = numel( xref );
    ny_ref = numel( yref );
    
    % Total length of a vectorized grid matrix
    %  (this is also the number of columns in M and rows in F)
    Nr = nx_ref*ny_ref;
    
    % Get cell width, note: cubic arrays must be uniform 
    dx = diff( xref ); 
    dy = diff( yref );
    
    % Tack on the last cell 
    % ( this isn't used b/c points using it should be culled)
    dx( end + 1 ) = dx( end );
    dy( end + 1 ) = dy( end );
    
    %% Reference Array Checks
    
    % Make sure points don't overlap
    if any( dx==0 ) || any( dy==0 )
        error('grid not distinct')
    end
    
    % Make sure grids are sorted
    if ~issorted( xref ) || ~issorted( yref )     
        error('grid is not sorted')
    end
    
    %% Global Indices & Local Coordinates
    
    % Get the Left-Hand Side index, points to 
    % where target resides in reference grid
    % (i.e. if xref(4:5)=[0.4,0.9] and xarray=0.7, returns ix = 4 )
    % [~,~,~, i_glbl_pntr, j_glbl_pntr ] = histcounts2( Xq, Yq, xref, yref );
    [~,~,i_glbl_pntr ] = histcounts( Xq, xref );
    [~,~,j_glbl_pntr ] = histcounts( Yq, yref );
    
    % Determine subset of indices that should not be used 
    %   ( histcounts return 0 for out of bounds)
    idx_null = ( i_glbl_pntr == 0 | j_glbl_pntr ==0 );
    
    % Give the bad cases an arbitrary valid LHS indices
    % (this preserves the target-array indices, should not 
    %  really affect performance if bad cases are rare)
    i_glbl_pntr( idx_null ) = 1;
    j_glbl_pntr( idx_null ) = 1;
    
    % Get the relative position in the [0,1x0,1] reference cell
    x(1,:) = ( Xq - xref( i_glbl_pntr ) ) ./ dx( i_glbl_pntr );
    y(1,:) = ( Yq - yref( j_glbl_pntr ) ) ./ dy( j_glbl_pntr );
    
    %% Cubic BC Handeling
    
    % Default use extrapolation
    bc_flag = false(4,1); % 1-4 are top, bot, left, right 
    
    % Build up cubic BC's
    if ( n_order==2 && ~isempty(cubic_bc_x) )
       
        % Determine which array is longer than nx/nyref
        if ( size(cubic_bc_x,1) == nx_ref+2 )
           nxref_bc = nx_ref+2;
           nyref_bc = ny_ref;
        elseif ( size(cubic_bc_y,1) == ny_ref+2 )
           nxref_bc = nx_ref;
           nyref_bc = ny_ref+2;
        else
            error('One Array of Cubic BCs must be n+2 (for corners)')
        end
        
        % Check inner sizes now 
        % BC notation: left x=0, right x=max, top y=max, bot y=0)
        % NOTE: bc_x contains i-top, j-top, i-bot, j-bot
        % NOTE: bc_y contains i-left, j-left, i-bot, j-bot
%         if ( size(cubic_bc_x,1)~= nxref_bc )
%            error('Incorrect inner dimension of cubic bc on top/bottom') 
%         elseif ( size(cubic_bc_y,1)~= nyref_bc )
%            error('Incorrect inner dimension of cubic bc on left/right') 
%         end
        
        % Create test of which walls to apply, turns true if not all zero
        bc_flag(1) = max( cubic_bc_x(:,1) ) > 0;
        bc_flag(2) = max( cubic_bc_x(:,3) ) > 0;
        bc_flag(3) = max( cubic_bc_y(:,1) ) > 0;
        bc_flag(4) = max( cubic_bc_y(:,3) ) > 0;        
        
        % check if any zero indicies in non-zero arrays
        if     ( min([cubic_bc_x(:,1);cubic_bc_x(:,2)])<1 && bc_flag(1) )
            error('Top BCs contain negative or zero indicies')
        elseif ( min([cubic_bc_x(:,3);cubic_bc_x(:,4)])<1 && bc_flag(2) )
            error('Bot BCs contain negative or zero indicies')
        elseif ( min([cubic_bc_y(:,1);cubic_bc_y(:,2)])<1 && bc_flag(3) )
            error('Left BCs contain negative or zero indicies')
        elseif ( min([cubic_bc_y(:,3);cubic_bc_y(:,4)])<1 && bc_flag(4) )
            error('Right BCs contain negative or zero indicies')
        end
        
    end
    
    %% Interpolation Coefficients
    % Need to get a square matrix C of coefficients that correspond to the
    % weights of the function points surrounding any given target x/y
    % this is stored in the inner index (rows) for each x/y target (column)
    %   Vectors A(x) and B(y) store the x/y-based weights that form C 
    
    % Linear (2x2 kernel) or nearest
    if ( n_order<2 )
        
        % if nearest, just round x and y
        if (n_order==0)
           x = round(x);
           y = round(y);
        end
        
        % Stencil size
        Ns = 4;
        
        % Weight Coefficients
        W(1,1,:) = (1.0-x).*(1.0-y);
        W(2,1,:) = (1.0-y).*x;
        W(1,2,:) = (1.0-x).*y;
        W(2,2,:) = x.*y;
        
        % Indices of the 2x2 kernel relative to the ix/iy global coordinate
        i_sten_mat(:,1) = [1:2,1:2]-1;
        j_sten_mat(:,1) = [1,1,2,2]-1;
        
    % Cubic (4x4 kernel)
    elseif (n_order==2)
        
        % Stencil size
        Ns = 16;
        
        % Substitution variables
        x2 = x.*x;
        y2 = y.*y;
        x2b = x-1.0;
        y2b = y-1.0;
        
        % Get A and B
        A(1,1,:) = -y2b.*y2b.*y;
        A(1,2,:) =  3.0.*y.*y2-5.0.*y2+2.0;
        A(1,3,:) = -3.0.*y.*y2+4.0.*y2+y;
        A(1,4,:) =  y2b.*y2;
        B(1,1,:) = -x2b.*x2b.*x;
        B(2,1,:) =  3.0.*x.*x2-5.0.*x2+2.0 ;
        B(3,1,:) = -3.0.*x.*x2+4.0.*x2+x;
        B(4,1,:) =  x2b.*x2;
        
        % Multiply by the common one-half factor
        A = A * 0.5; 
        B = B * 0.5;
        
        % Construct Matrix
        W = B.*A;
        
        % Apply 2nd order difference to extrapolate the outer edges of the
        % stencil when they are not available. This treatment ignores the
        % corners of the matrix. Deviations from interp2 are negligible,
        % less than <1% on average. 
        
        % Test and Apply Extrapolation Boundary condition on top wall
        top_idx(:,1) = find( j_glbl_pntr == Nr-1) ;
        if( ~bc_flag(1) )
            W(:,1,top_idx) = W(:,1,top_idx) + 0.5*W(:,4,top_idx);
            W(:,2,top_idx) = W(:,2,top_idx) - 2.0*W(:,4,top_idx);
            W(:,3,top_idx) = W(:,3,top_idx) + 2.5*W(:,4,top_idx);
            W(:,4,top_idx) = 0.0;
        end
        
        % Test and Apply Extrapolation Boundary condition on bottom wall
        bot_idx(:,1) = find( j_glbl_pntr == 1 );
        if( ~bc_flag(2) )
            W(:,2,bot_idx) = W(:,2,bot_idx) + 2.5*W(:,4,bot_idx);
            W(:,3,bot_idx) = W(:,3,bot_idx) - 2.0*W(:,4,bot_idx);
            W(:,4,bot_idx) = W(:,4,bot_idx) + 0.5*W(:,4,bot_idx);
            W(:,1,bot_idx) = 0.0;
        end
        
        % Test and Apply Extrapolation Boundary condition on left wall
        left_idx(:,1) = find( i_glbl_pntr == 1 );
        if( ~bc_flag(3) )
            W(2,:,left_idx) = W(2,:,left_idx) + 2.5*W(4,:,left_idx);
            W(3,:,left_idx) = W(3,:,left_idx) - 2.0*W(4,:,left_idx);
            W(4,:,left_idx) = W(4,:,left_idx) + 0.5*W(4,:,left_idx);
            W(1,:,left_idx) = 0.0;
        end
        
        % Test and Apply Extrapolation Boundary condition on right wall
        right_idx(:,1) = find( i_glbl_pntr == Nr-1 );
        if( ~bc_flag(4) )
            W(1,:,right_idx) = W(1,:,right_idx) + 0.5*W(4,:,right_idx);
            W(2,:,right_idx) = W(2,:,right_idx) - 2.0*W(4,:,right_idx);
            W(3,:,right_idx) = W(3,:,right_idx) + 2.5*W(4,:,right_idx);
            W(4,:,right_idx) = 0.0;
        end
        
        % Eliminate any boundaries where the BC is not needed because no
        % queries lie there
        bc_flag(1) = ~( bc_flag(1) && isempty( top_idx )   ); 
        bc_flag(2) = ~( bc_flag(2) && isempty( bot_idx )   ); 
        bc_flag(3) = ~( bc_flag(3) && isempty( left_idx )  ); 
        bc_flag(4) = ~( bc_flag(4) && isempty( right_idx ) ); 
        
        % Indices of the 4x4 kernel relative to the histcounts rules
        %  ( i.e. coordinates relative to i=1,j=1 being the LHS of x/y)
        i_sten_mat(:,1) = repmat(1:4,1,4)-2;
        j_sten_mat(:,1) = sort(repmat(1:4,1,4))-2;
        
    end
    
    % Now cull any C's that correspond to bad indices
    % (and replace with edgeval)
    W(:,:,idx_null) = 0.0; % sparse wont read zeros!
    
    % Reshape C to be Ns x Nq
    W = reshape(W, [Ns, Nq] );
    
    %% Global indices
       % Now get iloc(Ns,Nq) and jloc that contain the global
       % position of each stencil weight
    
    % Map the local stencil indices to global indices
    % glbl_mat i/j points to xref/yref indices for each query
    i_glbl_sten = i_glbl_pntr + i_sten_mat;
    j_glbl_sten = j_glbl_pntr + j_sten_mat;
    
    % Now fix boundary conditions for cubic case
    % NEED TO IMPROVE ANNOTATION
    if ( n_order==2 && any(bc_flag) )
        
        % Get X  modifier if corners are defined alonge the X axis
        % (top/bottom). If so, this shifts ij_bin to the correct indicies
%         px = 0;
%         if ( nxref_bc>nx_ref )
            px = 1;
            
%         end
        
        % Get Y modifier if corners are defined along the Y axis
        % (left/right) If so, this shifts ij_bin to the correct indicies
%         py = 0;
%         if ( nyref_bc>ny_ref )
            py = 1;
%         end           
            
        % Start with BC's defined along X (at Y=0 and Y=Ymax, or bottom and
        % top) use cubic_bc_x for bc 1 & 2, indexed by nx in the inner dim.
        % Assign columns of cubic_bc_x(:,1/3) and cubic_bc_x(:,2/4) to values
        % of i_glbl_sten and j_glbl_sten's first and last 4 stencil
        % coefficients. 
        
        % TOP CONDITION 
        % Use top_idx, which picks Y=Ymax cases out of i_glbl_pntr array
        % Use cubic_bc_x(:,1,2), 1 = X indicies to map to, 2 = Y indicies
        if( bc_flag(1) )
            
            % Location of top wall NY+1 node in every stencil with LHS = Ny-1
            idx_sten = 13:16; % Last 4 for Top-Wall, need to access V's outside
            
            % ij_bin advanced by 1 b/c user supplied 2 extra BC's on X for
            % the corners, meaning each index is shifted forwarded by 1
            ij_bin = i_glbl_sten( idx_sten, top_idx) + px; % i coordinate for all the j's at each bot_idx
            
            % Update the i (X) Bins
            i_glbl_sten( idx_sten, top_idx) = reshape( cubic_bc_x( ij_bin ,1), size(ij_bin));
            
            % Update the j (Y) bins
            j_glbl_sten( idx_sten, top_idx) = reshape( cubic_bc_x( ij_bin, 2), size(ij_bin));
            
        end
        
        % BOTTTOM CONDITION 
        % Use bot_idx, which picks Y=0 wall cases out of i_glbl_pntr array
        % Use cubic_bc_x(:,3,4), 3 = X indicies to map to, 4 = Y indicies
        if( bc_flag(2) )
            
            % Location of bottom-wall NY-1 node in every stencil with LHS = 1
            idx_sten = 1:4; % first 4 for bottom wall, need to access Y<0
            
            % IJ bin, positions along the cubic_bc_x reference array 
            ij_bin = i_glbl_sten( idx_sten, bot_idx) + px; % i coordinate for all the j's at each bot_idx
            
            % Update the i (X) Bins
            i_glbl_sten( idx_sten, bot_idx) = reshape( cubic_bc_x( ij_bin ,3), size(ij_bin));
            
            % Update the j (Y) bins
            j_glbl_sten( idx_sten, bot_idx) = reshape( cubic_bc_x( ij_bin, 4), size(ij_bin));
            
        end

        % Start with BC's defined along Y (at X=0 and X=Xmax, or left and
        % right. use cubic_bc_y for bc 3 & 4, indexed by ny in the inner
        % dim. Assign columns of cubic_bc_y(:,1/3) and cubic_bc_y(:,2/4) to
        % values of i_glbl_sten and j_glbl_sten respectivley,   stencil
        % coefficients.
        
        % LEFT CONDITION 
        % Use left_idx, which picks X=0 cases out of i_glbl_pntr array
        % Use cubic_bc_y(:,1/2), 1 = X indicies to map to, 2 = Y indicies
        if( bc_flag(3) )
            
            % Location of Left Wall Nx-1 in the local stencil
            idx_sten = [1,5,9,13]; 
            
            % ij_bin advanced by 1 b/c user supplied 2 extra BC's on X/Y for
            % the corners, meaning each index is shifted forwarded by 1
            ij_bin = j_glbl_sten( idx_sten, left_idx) + py; % i coordinate for all the j's at each bot_idx
            
            % Update the i (X) Bins
            i_glbl_sten( idx_sten, left_idx) = reshape( cubic_bc_y( ij_bin ,1), size(ij_bin));
            
            % Update the j (Y) bins
            j_glbl_sten( idx_sten, left_idx) = reshape( cubic_bc_y( ij_bin, 2), size(ij_bin));
            
        end
        
        % Right CONDITION 
        % Use right_idx, which picks X=xmax cases out of i_glbl_pntr array
        % Use cubic_bc_y(:,3/4), 3 = X indicies to map to, 4 = Y indicies
        if( bc_flag(4) )
            
            % Location of Right Wall Nx+1 in the local stencil
            idx_sten = [4,8,12,16]; 
            
            % ij_bin advanced by 1 b/c user supplied 2 extra BC's on X/Y for
            % the corners, meaning each index is shifted forwarded by 1
            ij_bin = j_glbl_sten( idx_sten, right_idx) + py; % i coordinate for all the j's at each bot_idx
            
            % Update the i (X) Bins
            i_glbl_sten( idx_sten, right_idx) = reshape( cubic_bc_y( ij_bin ,3), size(ij_bin));
            
            % Update the j (Y) bins
            j_glbl_sten( idx_sten, right_idx) = reshape( cubic_bc_y( ij_bin, 4), size(ij_bin));
            
        end
        

    elseif (n_order==2)
        

        
    end
    
    % Just change bad indicies, which should already be removed by
    % zeroing weights. If you dont have this, sub2ind fails
    i_glbl_sten( i_glbl_sten > nx_ref ) = nx_ref;
    i_glbl_sten( i_glbl_sten < 1      ) = 1;
    j_glbl_sten( j_glbl_sten > ny_ref ) = ny_ref;
    j_glbl_sten( j_glbl_sten < 1      ) = 1;
    
    %% Indexing the Sparse Matrix
       % We have iloc,jloc,W of size Ns x Nq,
       % corresponding to the global indices and weights
       % of the reference grid that are used
       % to construct any given x/y query. Now want to build the final
       % sparse matrix

    % Get row indices for each column of Nq
    I( 1:Ns, 1:Nq ) = repmat( 1:Nq, [Ns,1] ) ;
    
    % Get the column indices for each query (Ns,Nq)
    %  (based on the dimensions of the implied ny/nx grid)
    J = sub2ind( [nx_ref, ny_ref], i_glbl_sten, j_glbl_sten );
    
    % Create the sparse array
    M = sparse( I, J, W, Nq, Nr );
    
    %% Alternative Output
        
    % Define Output X and Y reference arrays, in the implied ndgrid form 
    % used internally
    if ( nargout>1 )
        
        % Copy X/Y to xout/yout
        xout = X;
        yout = Y;
        
        % If not already ndgrid, generate ndgrid
        if (mtype~=2)
            [xout,yout] = ndgrid(xref,yref);
        end
        
        % Then Flatten
        xout = xout(:);
        yout = yout(:);
        
    end
    
    % Define out-of-bounds pointer
    if ( nargout>3 )
        
        % Get inidices of out-of-bounds queries, vector form
        ob_pntr(:,1) = find(idx_null);
        
    end
    
end
