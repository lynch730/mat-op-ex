
close all
clear; clc;
addpath('core/');

% This script demonstrates the matrix extraction method of solving interp2
% Currently this is limited to only linear or cubic interpolation
% Currently both cases only accept regular Cartesian grids. 
% This case shows 7 test cases and their relative accuracy compared to
% interp2 with cubic interpolation
% NOTE: slight deviation on wall cases, due to uncertainty over interp2 BC.

%% Reference Data Setup

% Get reference data
nx = 18;
ny = 20;

% Define Grid Limits
xref(:,1) = linspace( -1.0, 1.0, nx );
yref(:,1) = linspace( -1.0, 1.2, ny ); 

% Create Mesh with Ngrid or meshgrid
[ X , Y ] = ndgrid(   xref, yref );
% [ X , Y ] = meshgrid(   xref, yref );

% Create Test Surface
Z = 3.0*peaks(X,Y);

%% Query Data Setup

% Interpolation Type (cubic or linear)
interp_type = 'cubic';

% Query Points

% Off-Grid Case
xquery(1,1) = xref(1)-0.05*range(xref);
yquery(1,1) = yref( round(ny/2.0) );

% Center
xquery(2,1) = 0.0;
yquery(2,1) = 0.0;

% Interior (normal) case
xquery(3,1) = xref(end) - 0.6*range(xref);
yquery(3,1) = yref(end) - 0.25*range(yref);

% Corner Case
xquery(4,1) = mean(xref(1:2));
yquery(4,1) = mean(yref(1:2));

% Wall Case
xquery(5,1) = xref(end) - 0.25*range(xref);
yquery(5,1) = mean(yref(1:2));

% On a line Case
xquery(6,1) = xref(end-3);
yquery(6,1) = mean(yref)-0.01*range(yref);

% On a Point Case
xquery(7,1) = xref(4);
yquery(7,1) = yref(end-4);


%% Getting Interpolations from matrix method

% Get size for reshaping after
nsize = size(xquery);

% Get Sparse Matrix of Remap Coefficients
tic
[M,xref_out,yref_out] = interp2_matrix( X, Y, xquery, ...
                                        yquery, interp_type );
toc

% Transpose if using meshgrid
mtype = mesh_type( X, Y );
if (mtype == 3)
   Z = Z'; 
end

% Get Values by matrix multiplication
Zval = M*Z(:);

% Reshape to original matrix
Zval = reshape(Zval,nsize);


%% Testing Against Interp2

% Recreate the meshgrid needed for classic interp2
[X2  ,Y2 ] = meshgrid( xref, yref );
Z2 = 3.0*peaks(X2,Y2);

% Run equivalent interp2
Ztest(:,1) = interp2( X2, Y2, Z2, xquery, yquery, interp_type);

% Percent Error
error_vals = abs((Zval-Ztest)./Ztest);


%% Plot Errors for each Example
 
% Open Figure
fig = figure(1);
fig.Units = 'inches';
fig.Position= [2,2,9,5];

% Array of Colors 
cc = [ 0.9047, 0.1918, 0.1988 ; ...
       0.2941, 0.5447, 0.7494 ; ...
       0.3718, 0.7176, 0.3612 ; ...
       1.0000, 0.5482, 0.1000 ; ...
       0.8650, 0.8110, 0.4330 ; ...
       0.6859, 0.4035, 0.2412 ; ...
       0.9718, 0.5553, 0.7741   ];

% Subplot for errors
subplot(1,4,1); hold on; grid on;

% Bar Chart
b=bar(100.*error_vals);
b.FaceColor = 'flat';

% Change colors to match scheme
for i = 1:length(xquery)
    b.CData(i,:) = cc(i,:);
end

% Plot Details
xlabel('Test Cases')
ylabel('Percent Error from Interp2')
title('Deviation')


%% Plot Stencils on the Z-reference surface

% Subplot ID
subplot(1,4,[2:4]); hold on; 

% Get location above surface to draw on
zmax = max(max(Z2))*1.05;

% Plot Surface z
pp(1) = surf(X2,Y2,Z2);
pp(1).FaceAlpha = 0.8;
view(0,90);colorbar;
colormap gray;

% Plot Each Stencil in example list
for i = 1:length(xquery)
    
    % Plot Point
    pp(i+1)= plot3(xquery(i),yquery(i),zmax,'x',...
          'MarkerSize',12 , ...
          'LineWidth' ,2.5, 'Color', cc(i,:) );
    
    % Plot Stencil Coordinates, if they are valid
    if ~isnan( Zval(i) )
        
        % Get index of column in smat to pull from for each query
        idx = find(M(i,:));
        
        % Arbitrary Z function to make surface visible above mesh
        zvector = repmat( zmax, length(idx) ,1 );
        
        % plot stencil points as scatter
        plot3(xref_out(idx),yref_out(idx),zvector, ...
              '.','MarkerSize',20,...
              'Color',cc(i,:)  );
          
    end
    
end

% Legend for both subplots
legend( pp(1:end), 'Z Surface',...
        '1) Off-Grid',...
        '2) Center',...
        '3) Interior', ...
        '4) Corner', ...
        '5) Near-Wall', ...
        '6) On Line', ...
        '7) On Node',...
        'NumColumns',1 )
    
% plot Details
xlim( [ min( [xref;xquery] ), max( [xref ; xquery] ) ] )
ylim( [ min( [yref;yquery] ), max( [yref;yquery  ] ) ] )
xlabel('X')
title('Stencils')
set(findobj(gcf,'type','axes'),'FontName','Arial',...
    'FontSize',10,'FontWeight','Bold');

