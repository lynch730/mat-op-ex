function [ M, ob_pntr ] = interp1_matrix( xref, Xq, varargin )
%% INTERP1_MATRIX, RECOVERING THE MATRIX FROM INTERP1
%     
%     M = interp1_matrix(xref,Xq) returns the underlying matrix M
%     representing the interp1 operator for a known xref grid and set of a
%     set of query points Xq. Matrix M is of size M(Nq,Nr), with Nq rows
%     of vectorized queries Xq(:), and Nr columns of vectorized
%     reference grid points xref(:).
%     
%     M can subsequently be used to obtain a vector P(Nq,1) of
%     Nq interpolation queries for a reference function F over xref:
%        P( 1:Nq,1 ) = M( 1:Nq, 1:Nr ) * F( 1:Nr,1 )
%     
%     While M could in theory be constructed to preserve the shapes of xref
%     and Xq in a high-order tensor that is contracted over
%     F(xref), interp1_matrix vectorizes to maintain compatibility with
%     MATLAB matrix multiplication. It is left to the user to reshape M if
%     that is desired.
% 
%
%% FUNCTION OPTIONS
% 
%     M = interp1_matrix(xref,Xq,METHOD) specifies interpolation method,
%     which is linear by default. Other options are nearest and cubic.
%     Results should be identical to interp1, with the exception of cubic
%     cases where the stencil requires cells outside the grid.
%     interp1_matrix uses a symmetry condition for those cells, which
%     incurs deviations in these cells of a few percent compared to
%     interp1.
%     
%       'nearest' - nearest neighbor interpolation 
%       'linear'  - linear interpolation 
%       'cubic'   - cubic convolution interpolation
%
%    [ M, OB_PNTR ] = interp1_matrix(...) returns an array of
%    the M matrix's row indicies corresponding to query points that are out
%    of bounds. These rows are left sparse. 
%    
%% SOLUTION TECHNIQUE & SOURCE
%    
%    This solution is based heavily on the suggested technique in
%     <https://www.mathworks.com/matlabcentral/answers/573703-get-interpolation-transfer-relation-matrix-instead-of-interpolated-values>
%    which is expanded upon to include the cubic case
%    
%    The stencil "sten" is the grid points associated with individual interpolations.
%    Global "glbl" refers to indices/coordinates in the frame of the larger
%    xref reference grid "ref".
%    
%    The general approach for each query is:
%       (1) Build the arrays i_glbl_pntr (integer arrays which
%           point each query to the nearest left-hand-side grid point in
%           xref)
%       (2) Find the appropriate stencil (2x1 for linear and nearest,
%           4x1 cubic) of weights W, based on cell-relative position "x"
%       (3) Compute i_glbl_sten which map each W to the reference grid.
%       (5) Construct the sparse matrix from W,I,J, where I is the row
%           index of M, the same for all stencil weights in each query.
%    The above (1-5) are all simultaneously performed for each element in
%    Xq, indexed on the outside of all of the above operations.
%    
%    Though not wildly slow, there is likely room to optimize the current
%    technique. However user accessibility is perhaps more important here.
%    The primary utility in obtaining the matrix is for minimizing the
%    overhead associated with repeated interp1 calls on a fixed grid, where
%    the M matrix need only be computed once. Note, this can be done to
%    some extent for conventional CPU operations by creating
%    P=griddedInterpolant (the underlying code interp1 wraps around), and
%    updating only the F.Values to avoid recomputing the weights.
%    Critically, this approach currently cannot be replicated for
%    gpuArrays, where such optimizations would be most useful. Of course,
%    interp2_matrix is still useful for analysis of Jacobians and other
%    applications that griddedInterpolant cannot reproduce.
%    

    %% Process Optional Arguments
    
    % Specify input type, if provided
    method = 'linear';
    if ( nargin > 2 )
       method = varargin{1};
    end
    
    %% Check the Interpolation Method
    
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
    
    %% Pre-process the Query Data
    
    % Reshape to linear, with queries stored along outer index (columns)
    Xq = reshape( Xq, 1, [] );
    
    % Get the length of the linear arrays
    Nq = numel( Xq );
    
    %% Pre-process the Reference Data
    
    % Check if successful
    if (~isvector(xref))
        error('xref is not the correct format')
    end
    
    % Reshape to vector, stored like target arrays (outer, in columns)
    xref = reshape( xref, 1, [] );
    
    % Get lengths of the reference arrays, these can be different
    Nr = numel( xref );
    
    % Get cell width, note: cubic arrays must be uniform 
    dx = diff( xref ); 
    
    % Tack on the last cell 
    % ( this isn't used b/c points using it should be culled)
    dx( end + 1 ) = dx( end );
    
    %% Reference Array Checks
    
    % Make sure points don't overlap
    if any( dx==0 )
        error('grid not distinct')
    end
    
    % Make sure grids are sorted
    if ~issorted( xref ) 
        error('grid is not sorted')
    end
    
    %% Global Indices & Local Coordinates
    
    % Get the Left-Hand Side index, points to 
    % where target resides in reference grid
    [~,~,i_glbl_pntr ] = histcounts( Xq, xref );
    
    % Determine subset of indices that should not be used 
    %   ( histcounts return 0 for out of bounds)
    idx_null = ( i_glbl_pntr == 0 );
    
    % Give the bad cases an arbitrary valid LHS indices
    % (this preserves the target-array indices, should not 
    %  really affect performance if bad cases are rare)
    i_glbl_pntr( idx_null ) = 1;
    
    % Get the relative position in the [0,1x0,1] reference cell
    x(1,:) = ( Xq - xref( i_glbl_pntr ) ) ./ dx( i_glbl_pntr );
    
    %% Interpolation Coefficients
    % Need to get a square matrix C of coefficients that correspond to the
    % weights of the function points surrounding any given target x
    % this is stored in the inner index (rows) for each x target (column)
    
    % Linear (2x1 kernel) or nearest
    if ( n_order<2 )
        
        % if nearest, just round x and y
        if (n_order==0)
           x = round(x);
        end
        
        % Stencil size
        Ns = 2;
        
        % Weight Coefficients
        W(1,:) = 1.0-x;
        W(2,:) = x;
        
        % Indices of the 2x1 kernel relative to the ix global coordinate
        i_sten_mat(:,1) = [0,1];
        
    % Cubic (4x4 kernel)
    elseif (n_order==2)
        
        % Stencil size
        Ns = 4;
        
        % X temporary
        x2 = x.*x;
        x2b = x-1.0;
        
        % Get A and B
        W(1,:) = -x2b.*x2b.*x;
        W(2,:) =  3.0.*x.*x2-5.0.*x2+2.0 ;
        W(3,:) = -3.0.*x.*x2+4.0.*x2+x;
        W(4,:) =  x2b.*x2;
        
        % Multiply by the common one-half factor
        W = W * 0.5;        
        
        % Apply Boundary condition on left wall
        wall = (i_glbl_pntr == 1);
        W(2,wall) = W(2,wall) + 2.5*W(4,wall);
        W(3,wall) = W(3,wall) - 2.0*W(4,wall);
        W(4,wall) = W(4,wall) + 0.5*W(4,wall);
        W(1,wall) = 0.0;
        
        % Apply Boundary condition on right wall
        wall = (i_glbl_pntr == Nr-1);
        W(1,wall) = W(1,wall) + 0.5*W(4,wall);
        W(2,wall) = W(2,wall) - 2.0*W(4,wall);
        W(3,wall) = W(3,wall) + 2.5*W(4,wall);
        W(4,wall) = 0.0;
        
        % Indices of the 4x1 kernel relative to the histcounts rules
        %  ( i.e. coordinates relative to i=1 being the LHS of x/y)
        i_sten_mat(:,1) = [-1,0,1,2];
        
    end
    
    
    %% Global indices
       % Now get iloc(Ns,Nq) and jloc that contain the global
       % position of each stencil weight
    
    % Map the local stencil indices to global indices
    % glbl_mat i points to xref indices for each query
    i_glbl_sten = i_glbl_pntr + i_sten_mat;
    
    % Now cull any C's that correspond to bad indices
    % (and replace with edgeval)
    W(:,idx_null) = 0.0;
    
    % Apply a symmetry condition to reference points that 
    % lie outside the domain. This assumes zero slope at domain edge
    %   ( In practice, just need to inc/decrement any references to 0/end+1)
    %   ( Only do this on CUBIC type, stencil for linear doesn't ever need end+1)
    if ( n_order==2 )
        i_glbl_sten( i_glbl_sten > Nr ) = Nr;
        i_glbl_sten( i_glbl_sten < 1  ) = 1;
    end
    
    %% Indexing the Sparse Matrix
       % We have iloc of size Ns x Nq,
       % corresponding to the global indices and weights
       % of the reference grid that are used
       % to construct any given x query. Now want to build the final
       % sparse matrix

    % Get row indices for each column of Nq
    I( 1:Ns, 1:Nq ) = repmat( 1:Nq, [Ns,1] ) ;
    
    % Create the sparse array
    M = sparse( I, i_glbl_sten, W, Nq, Nr );
    
    %% Alternative Output
        
    % Return pointer array of indicies in rows of M that were nullfied due
    % to out-of-bounds queries
    if ( nargout>1 )
        ob_pntr(:,1) = find(idx_null);
    end
    
    
end

