
close all
clear; clc;
if (ispc)
   error('If using windows, add path to source here')
%   addpath('...\source\')
else
   addpath('../../source/')
end

% This script demonstrates the matrix extraction method on an interp2 case
% where the entire surface is shifted in the X&Y-axis by 1/10 the X/Y
% range. This requires remapping every cell, and handling edge cases. Spy
% analysis of the density of the sparse matrix is included
%
% NOTE: regions of shift plot contain out-of-bound NaNs, and are not
% plotted

%% Reference Data Setup

% Get reference data
nx = 18;
ny = 20;

% Define Grid Limits
xref(:,1) = linspace( -1.0, 1.0, nx );
yref(:,1) = linspace( -1.0, 1.2, ny ); 

% Create Mesh with Ngrid or meshgrid
[ X , Y ] = ndgrid(   xref, yref );

% Create Test Surface
Z = 3.0*peaks(X,Y);


%% Problem Setup

% Interpolation Type (cubic or linear)
interp_type = 'cubic';

% Setup queries to seek out the point down and to the left
xquery = X - 0.1*range(xref);
yquery = Y - 0.1*range(yref);


%% Getting Interpolations from matrix method

% Get size for reshaping after
nsize = size(xquery);

% Get Sparse Matrix of Remap Coefficients
[ M, ~, ~, ob_pntr ] = interp2_matrix( X, Y, xquery, yquery, interp_type );

% Get Values by matrix multiplication
Z2 = M*Z(:);

% Reshape to original matrix
Z2 = reshape(Z2,nsize);

% Apply Constant Condition on Out-of-bounds queries
Z2( ob_pntr ) = mean(mean(Z)); 

%% Plot Before and After
 
% Open Figure
fig = figure(1);
fig.Units = 'inches';
fig.Position= [2,2,9,4];

% Original
subplot(1,2,1); 
pp(1) = surf(X,Y,Z);
view(0,90);colorbar;
title('Original')
xlabel('X')
ylabel('Y')
caxis([ min(min(Z)) , max(max(Z)) ])
xlim(  [ min( xref ), max( xref )  ] )
ylim(  [ min( yref ), max( yref )  ] )

% Shifted
subplot(1,2,2); 
pp(2) = surf(X,Y,Z2);
view(0,90);colorbar;
colormap gray;
title('Shifted')
xlabel('X')
ylabel('Y')
caxis( [ min(min(Z)) , max(max(Z)) ] )
xlim(  [ min( xref ), max( xref )  ] )
ylim(  [ min( yref ), max( yref )  ] )
set(findobj(gcf,'type','axes'),'FontName','Arial',...
     'FontSize',10,'FontWeight','Bold');
 
 
%% Plot Spy
 
% Open Figure
fig = figure(2); grid on;
fig.Units = 'inches';
fig.Position= [2,2,5,5];

% plot
spy(M)
title('M Interpolation Matrix')
xlabel('N_x x N_y Reference Values')
ylabel('N_x x N_y Query Values')

