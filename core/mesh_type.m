%﻿% Auxiliary function, determine format of X and Y matrices
function [mtype, xarray, yarray ]= mesh_type( X, Y )
    
    % Returns mtype = 
    %  0, not a valid mesh type
    %  1, both X and Y are vectors
    %  2, both X and Y are ndgrid format
    %  3, both X and Y are meshgrid format
    
    % Default mtype to fail
    mtype = 0;
    xarray = NaN;
    yarray = NaN;
    
    % Test if both are the same size
    if  size( X ) ~= size( Y ) 
        return
    end
    
    % Test if vectors
    if isvector( X )
        mtype = 1;
        
    % test if matrices
    elseif ismatrix( X ) 
        
        % Tolerance for equality
        tol = 1e-10;

        % Determine which dimensions of X and Y repeat (logical)
        x_row_repeat = abs( ( X(end,1)-X(1,1) ) ) < tol; 
        x_col_repeat = abs( ( X(1,end)-X(1,1) ) ) < tol; 
        y_row_repeat = abs( ( Y(end,1)-Y(1,1) ) ) < tol; 
        y_col_repeat = abs( ( Y(1,end)-Y(1,1) ) ) < tol; 
        
        % Test if in ndgrid format 
        %   ( X-varies along dim=1 (rows), Y-varies along dim=2 (col))
        if ( x_col_repeat && y_row_repeat && ... 
            ~x_row_repeat && ~y_col_repeat )    

            % both X and Y are ndgrid forms
            mtype = 2;

        % Test if in meshgrid format 
        %   ( X-varies along dim=2 (col), Y-varies along dim=1 (row))
        elseif ( x_row_repeat && y_col_repeat && ... 
                ~x_col_repeat && ~y_row_repeat )   

            % both X and Y are meshgrid forms
            mtype = 3;

        end
        
    end
    
    %% Optional creation of the 1D grid vectors
    if nargout>1
        if (mtype == 1)  % vectors, copy directly
            xarray = X;
            yarray = Y;
        elseif (mtype == 2)  % ndgrid format
            xarray = X(:,1); % each row in 1st column 
            yarray = Y(1,:); % each column in 1st row
        elseif (mtype == 3)  % meshgrid format
            xarray = X(1,:); % each row in 1st column 
            yarray = Y(:,1); % each column in 1st row
        end
    end
    
end
