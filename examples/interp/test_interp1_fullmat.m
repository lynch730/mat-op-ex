
close all
clear; clc;
if (ispc)
   error('If using windows, add path to source here')
%   addpath('...\source\')
else
   addpath('../../source/')
end

% This script demonstrates the matrix extraction method on an interp1 case
% where the entire line is shifted in the X-axis by 1/10 the X range.
% This requires remapping every cell, and handling edge cases. Spy
% analysis of the density of the sparse matrix is included
%
% NOTE: regions of shift plot contain out-of-bound NaNs, and are not
% plotted

%% Reference Data Setup

% Get reference data
nx = 18;

% Define Grid Limits
xref(:,1) = linspace( -2.0, 3.0, nx );

% Create Test Surface
Z(:,1) = 10.0.*xref.*xref + 0.8.*xref + 1.2;


%% Problem Setup

% Interpolation Type (cubic or linear)
interp_type = 'cubic';

% Setup queries to seek out the point down and to the left
xquery = xref - 0.1*range(xref);


%% Getting Interpolations from matrix method

% Get size for reshaping after
nsize = size(xquery);

% Get Sparse Matrix of Remap Coefficients
[M, ob_pntr] = interp1_matrix( xref, xquery, interp_type , 35.0);

% Get Values by matrix multiplication
Zval = M*Z;

% Apply a condition for out-of-bounds query
Zval( ob_pntr ) = Zval(  find( (Zval>0.0), 1 ) ); 

%% Plot Before and After
 
% Open Figure
fig = figure(1);
fig.Units = 'inches';
fig.Position= [2,2,9,4];

% Original
subplot(1,2,1); grid on; hold on;
pp(1) = plot(xref,Z,'k.-');
title('Original')
xlabel('X')
ylabel('F')
xlim(  [ min( xref ), max( xref )  ] )
ylim(  [ min( Z    ), max( Z    )  ] )

% Shifted
subplot(1,2,2); grid on; hold on;
pp(2) = plot(xref,Zval,'k.-');
title('Shifted')
xlabel('X')
ylabel('F')
xlim(  [ min( xref ), max( xref )  ] )
ylim(  [ min( Z    ), max( Z    )  ] )

% Bold
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
xlabel('Nr Reference Values')
ylabel('Nq Query Values')

