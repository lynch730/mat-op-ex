
close all
clear; clc;
if (ispc)
   error('If using windows, add path to source here')
%   addpath('...\source\')
else
   addpath('../../source/')
end

% This script demonstrates the matrix extraction method of solving interp1
% Currently this is limited to only linear or cubic interpolation
% Currently both cases only accept regular cartesian grids. 
% This case shows 7 test cases and their relative accuracy compared to
% interp2 wtih cubic interpolation
% NOTE: slight deviation on wall cases, due to uncertainty over interp2 BC.

%% Reference Data Setup

% Set number of grid points
nx = 25;

% Define Grid Limits
xref(:,1) = linspace( -2.0, 3.0, nx );

% Create Test Surface
Z(:,1) = 10.0.*xref.*xref + 0.8.*xref + 1.2;

%% Query Data Setup

% Interpolation Type (cubic or linear)
interp_type = 'cubic';

% Query Points

% Off-Grid Case
xquery(1,1) = xref(1)-0.05*range(xref);

% Interior (normal) case
xquery(2,1) = xref(end) - 0.2*range(xref);

% % Corner Case
xquery(3,1) = mean(xref(end-1:end));
% 
% % On a point Case
xquery(4,1) = xref(4);

%% Getting Interpolations from matrix method

% matrix operator
M = interp1_matrix( xref, xquery, interp_type );

% interpolation by matrix multiplication
Zval = M*Z;


%% Testing Against Interp1

% Run equivalanet interp2
Ztest(:,1) = interp1( xref, Z, xquery, interp_type);

% Percent Error
error_vals = abs((Zval-Ztest)./Ztest);

%% Plot Errors for each Example
 
% Open Figure
fig = figure(1); 
fig.Units = 'inches';
fig.Position= [2,2,9,5];

% Array of Colors (from linspecer)
cc = [ 0.3467, 0.5360, 0.6907; ...
       0.9153, 0.2816, 0.2878; ...
       0.4416, 0.7490, 0.4322; ...
       1.0000, 0.5984, 0.2000  ];

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
ylabel('Percent Error from Interp1')
title('Deviation')


%% Plot Stencils on the Z-reference surface

% Subplot ID
subplot(1,4,2:4); grid on; hold on;

% Plot Line
pp(1) = plot(xref,Z,'k.-');

% Plot Each Stencil in example list
for i = 1:length(xquery)
    
    % Plot Point
    pp(i+1)= plot( xquery(i),Zval(i),'x',...
          'MarkerSize',12 , ...
          'LineWidth' ,2.5, 'Color', cc(i,:) );
    
    % Plot Stencil Coordinates, if they are valid
    if ~isnan( Zval(i) )
        
        % Get index of column in smat to pull from for each query
        idx = find(M(i,:));
        
        % plot stencil points as scatter
        plot( xref(idx),Z(idx), '.','MarkerSize',20,...
              'Color',cc(i,:)  );
          
    end
    
end

% Legend for both subplots
legend( pp(1:end), 'Z Solution',...
        '1) Off-Grid',...
        '2) Interior',...
        '3) Corner', ...
        '4) On Point',...
        'Location','NorthWest')
    
% plot Details
xlim( [ min( [xref;xquery] ), max( [xref ; xquery] ) ] )
xlabel('X')
title('Stencils')
set(findobj(gcf,'type','axes'),'FontName','Arial',...
    'FontSize',10,'FontWeight','Bold');

