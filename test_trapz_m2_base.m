
close all
clear; clc;
addpath('core/');

% This script demonstrates the matrix extraction method of solving trapz
% using matrix methods, applied to a problem of 2D grids. 

%% Reference Data Setup

% Get reference data
nx = 22;
ny = 19;

% Define Grid Limits
xref(:,1) = linspace( 0.0, 3.0, nx );
yref(:,1) = linspace( -3.0, -2.0, ny ); 

% Create Mesh with Ngrid or meshgrid
[ X , Y ] = ndgrid(   xref, yref );

% Create Test Surface
Z = -1.0 + 3.0*peaks(X,Y);

%% Test X Integration

% Dimension 
dim = 1;

% Trapz Test with Matrix Method
[Mx] = trapz_matrix2( X , Y,  dim, 'trapz' );
Zval_x = Mx*(Z(:));

% Composite Simpsons Rule with Matrix Method
[Mx] = trapz_matrix2( X , Y,  dim, 'simpsons' );
Zval_x_simp = Mx*(Z(:));

% Matlab Test
Ztest_x(:,1) = trapz( xref, Z, dim);

% Matlab Fine Grid
xref_f(:,1) = linspace( 0.0, 3.0, 1000 );
yref_f(:,1) = linspace( -3.0, -2.0, 1000 ); 
[ Xf , Yf ] = ndgrid(   xref_f, yref_f );
Zf = -1.0 + 3.0*peaks(Xf,Yf);
zt_x_f(:,1) = trapz( xref_f, Zf, dim);
Ztest_x_f(:,1) = interp1( yref_f, zt_x_f, yref  );

%% Test Y Integration

% Dimension 
dim = 2;

% Get Sparse Matrix of Remap Coefficients
[My] = trapz_matrix2( X , Y,  dim, 'trapz' );
Zval_y = My*(Z(:));

% Get Sparse Matrix of Remap Coefficients
[My] = trapz_matrix2( X , Y,  dim, 'simpsons' );
Zval_y_simp = My*(Z(:));

% Matlab Test
Ztest_y(:,1) = trapz( yref, Z, dim);

% Matlab Fine Grid
zt_y_f(:,1) = trapz( yref_f, Zf, dim);
Ztest_y_f(:,1) = interp1( xref_f, zt_y_f, xref  );

%% Errors

% % Percent Error
error_ztest_x     = abs( ( Ztest_x     - Ztest_x_f ) ./ Ztest_x_f );
error_zval_x      = abs( ( Zval_x      - Ztest_x_f ) ./ Ztest_x_f );
error_zval_x_simp = abs( ( Zval_x_simp - Ztest_x_f ) ./ Ztest_x_f );

error_ztest_y     = abs( ( Ztest_y     - Ztest_y_f ) ./ Ztest_y_f );
error_zval_y      = abs( ( Zval_y      - Ztest_y_f ) ./ Ztest_y_f );
error_zval_y_simp = abs( ( Zval_y_simp - Ztest_y_f ) ./ Ztest_y_f );

%% Plot Errors 
 
% % Open Figure
fig = figure(1); 
fig.Units = 'inches';
fig.Position= [2,2,9,5];

% Get min/max error in both frames
yall =[ error_zval_x_simp ; error_zval_x ; error_ztest_x ; ...
        error_zval_y ; error_ztest_y;  error_zval_y_simp ] ;
ymin = 100.0*min(yall);
ymax = 100.0*max(yall);

% Plot Dim1
subplot(1,2,1); 
grid on; hold on;
plot(yref,100.0.*error_zval_x,'-r','LineWidth',1.5)
plot(yref,100.0.*error_zval_x_simp,'-b','LineWidth',1.5)
plot(yref,100.0.*error_ztest_x,'xk','LineWidth',1.5)
xlabel('Y-Axis')
ylabel({'Percent Error','(Relative to trapz w/n=1000)'})
title('X-Integral Errors (DIM=1)')
ylim([ymin ymax])
set(gca,'YScale','log')

subplot(1,2,2); 
grid on; hold on;
plot(xref,100.0.*error_zval_y,'-r','LineWidth',1.5)
plot(xref,100.0.*error_zval_y_simp,'-b','LineWidth',1.5)
plot(xref,100.0.*error_ztest_y,'xk','LineWidth',1.5)
xlabel('X-Axis')
title('Y-Integral Errors (DIM=2)')
ylim([ymin ymax])
set(gca,'YScale','log')

% Legend
legend('Matrix2 - Trapezoidal','Matrix2 - Simpsons','MATLAB - trapz',...
      'Location','south' )

set(findobj(gcf,'type','axes'),'FontName','Arial',...
    'FontSize',10,'FontWeight','Bold');


%% Plot Integrals

% Open Figure
fig2 = figure(2); 
fig2.Units = 'inches';
fig2.Position= [2,2,9,5];

% Marker size
msize = 8;

% Plot Dim1
subplot(1,2,1); grid on; hold on;
plot(yref_f,zt_x_f,'-k','LineWidth',1.5)
plot(yref,Ztest_x,'og','MarkerSize',1.5*msize,'LineWidth',1.5)
plot(yref,Zval_x,'rx','MarkerSize',1.5*msize,'LineWidth',1.5)
plot(yref,Zval_x_simp,'b+','MarkerSize',msize,'LineWidth',1.5)
xlabel('Y-Axis')
ylabel('Integral')
title('X-Integral (DIM=1)')

% Plot Dim2
subplot(1,2,2); grid on; hold on;
plot(xref_f,zt_y_f,'-k','LineWidth',1.5)
plot(xref,Ztest_y,'og','MarkerSize',1.5*msize,'LineWidth',1.5)
plot(xref,Zval_y,'rx','MarkerSize',1.5*msize,'LineWidth',1.5)
plot(xref,Zval_y_simp,'b+','MarkerSize',msize,'LineWidth',1.5)
xlabel('X-Axis')
title('Y-Integral (DIM=2)')
legend('MATLAB - trapz (N=1000)','MATLAB - trapz',...
       'Matrix2 - Trapezoidal','Matrix2 - Simpsons',...
       'Location','southeast' )

% Set Bold
set(findobj(gcf,'type','axes'),'FontName','Arial',...
    'FontSize',10,'FontWeight','Bold');

%% Plot Spy

% Create figure
fig2 = figure(3); 
fig2.Units = 'inches';
fig2.Position= [2,2,9,5];

% Plot X-Integral
subplot(2,1,1)
spy(Mx)
title('X-Integral (DIM=1)')
ylabel('Y Result (Ny)')

% Plot Y-Integral
subplot(2,1,2)
spy(My)
title('Y-Integral (DIM=2)')
ylabel('X Result (Nx)')
xlabel('Number of X/Y Points (Nx*NY)')

% Set Bold
set(findobj(gcf,'type','axes'),'FontName','Arial',...
    'FontSize',10,'FontWeight','Bold');

