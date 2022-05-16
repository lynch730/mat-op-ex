function [ M ] = trapz_matrix1(  X, options )
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
    
    arguments
        X (1,:) double {mustBeVector(X), mustBeNumeric, issorted(X), mustBeReal} 
        options.M (:,:) double
        options.method (1,:) char {mustBeMember(...
                            options.method,{'simpsons','trapz'})} = 'trapz'
        options.interval (:,2) {mustBeInteger(options.interval)}
    end
    
%% Pre-process the Reference Data
    
    % Reshape to vector, stored like target arrays (outer, in columns)
    xref = reshape( X, [], 1 );
    
    % Get lengths of the reference arrays, these can be different
    Ncol = numel( xref );
    
    % Get cell width, note: cubic arrays must be uniform 
    dx(1,:) = diff( xref ); 
    
    % Tack on the last cell
    % ( this isn't used b/c points using it should be culled)
    dx( end + 1 ) = dx( end );
    
    % Make sure points don't overlap
    if any( dx==0 )
        error('grid not distinct')
    end
    
%     % Make sure grids are sorted
%     if ~issorted( xref )   
%         error('grid is not sorted')
%     end
    

   %%  Create Integration Coefficients V 

   % Check M input
   if isfield(options, 'M')
       if size(options.M, 2) ~= size(X, 2)
           error('Columns of M do not match X')           
       end
        M = options.M;
   else
       M = ones(1, Ncol);
   end
   Nrow = size(M, 1);

   % Check 
   if isfield(options, 'interval')
       if size(options.interval, 1) > 1 && Nrow==1
           Nrow = size(options.interval, 1);
       elseif size(options.interval, 1) ~= Nrow
           error('Rows of Interval do not match M')           
       end
       istart = options.interval(:, 1);
       iend = options.interval(:, 2);
       if max(iend)>Ncol || min(istart)<1
            error('interval contains OOB points')
       end
   else
        istart = ones(Nrow, 1);
        iend = istart + Ncol - 1;
   end

    % Multiply by axis delta
    M = M .* dx;
    
    % NOTE:
    % 5/13/22 - Added ad hoc rule that NaN's in M indicate regions where 
    % integration is not to occur, such that we "reset" integration bounds
    % to points adjacent. Implemented in trapz1 only, for the purpose
    % of allowing inerp+trap operations over restrctied finite volumes

    % Fork Integration Method
    if options.method == "trapz"
        
        % TrapZ Coefficients
        sz = [Nrow, Ncol];
        NrowM = size(M,1);
        Ncol_est = 10; % estimated number of columns used in each row
        i = zeros(Nrow*Ncol_est,1);
        j = zeros(Nrow*Ncol_est,1);
        v = zeros(Nrow*Ncol_est,1);
        cnt = 1;
        for irow =1:Nrow
            jnew = istart(irow):iend(irow);
            N = numel(jnew);
            Mvec = M(min(NrowM, irow), jnew);
            Mvec([1,N]) = Mvec([1,N]) *0.5;
            cnt2 = cnt+N-1;
            i(cnt:cnt2) = irow;
            j(cnt:cnt2) = jnew;
            v(cnt:cnt2) = Mvec;
            cnt = cnt2 + 1;
        end
        i(cnt:end) = [];
        j(cnt:end) = [];
        v(cnt:end) = [];
        M = sparse(i, j, v, sz(1), sz(2)); 
   
    elseif options.method == "simpsons"
        
        % TEMPORARY, FIND A VECOTRIZED FORM FOR THIS LOOP
        for i = 1:Nrow

            % Indicies to use
            iend_tmp = iend(i);
            ind = istart(i):iend_tmp;
            M(setdiff(1:Ncol, ind)) = 0.0;

            % If even grid points (odd panels)
            modified = false;
            if (mod(numel(ind), 2) == 0)
                modified = true;
                assert( numel(ind) > 4, 'Grid X too small (<4) for simpsons method')
                 
                % Reset limits for Simpson
                Nshift = 3;
                 iend_tmp = iend_tmp - Nshift;
                 ind = istart(i):iend_tmp;

                 % Save state of the last 4 terms
                 Vend(1, 1:Nshift+1) = M(i, iend_tmp+1:iend(i));
                 
            end
            
            % Simpsons Common Factor
            M(i,ind) = M(i,ind)./3.0;
            
            % Simpsons Multipliers
            M( i, ind(2:2:end-1)  ) = 4.0 .* M( i, ind(2:2:end-1) ) ;
            M( i, ind(3:2:end-2)  ) = 2.0 .* M( i, ind(3:2:end-2) ) ;
            
            % Finish Even case with Simpsons 3/8
            if modified
                
                % Overlap term used in both stencils
                M(i, ind(end))   = M(i, ind(end))   + 0.375*Vend(1, 1);
                
                % Remaining terms 
                M(i, ind(end)+1) = 1.1250*Vend(1, 2);
                M(i, ind(end)+2) = 1.1250*Vend(1, 3);
                M(i, ind(end)+3) =  0.375*Vend(1, 4);
                
            end

        end
        
    end
    
    % Return Sparse
    M = sparse(M);
    
end

