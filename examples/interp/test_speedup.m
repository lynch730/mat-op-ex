
close all
clear; clc;
if (ispc)
   error('If using windows, add path to source here')
%   addpath('...\source\')
else
   addpath('../../source/')
end

% This script tests the speedup of interp2_matrix for 2D matrix
% interpolation on a fixed grid, relative to interp2, griddedInterpolant,
% and with a GPU, if available. The script tests a range of grid sizes,
% meant to show the relative CPU/GPU scaling. This will take a long time,
% depending on the hardware and max grid size. Single GPU memory limits may
% require chaging max grid size

%% Hardcode CPU/GPU info
cpu_name = 'Intel i7-5820K';
gpu_name = 'NVIDIA GTX 1070';

% Number of iterations to time ( more will have smoother data )
% More iteraetions will favor schemes that don't remesh. 
ni = 1000;
min_grid = 10;
max_grid = 600;

% Loop number of grid points trials (log-spaced)
ns = 20;
nsa = linspace( log10(min_grid), log10(max_grid) , ns );
nsa = round(10.0.^nsa);

% Create Timer
timer = zeros(6,ns); 
timer_relative = timer;

% Open device, if possible
gpu_flag = false;
if (parallel.gpu.GPUDevice.isAvailable)
    
    % Open Device
    gpu_flag = true;
    gpu_device = gpuDevice;
    reset(gpu_device); 
    
end

% Main Loop
for itrial = 1:ns
    
    % Set New Ngrid
    nx = nsa(itrial);
    ny = nx;

    %% Setup Reference Data 

    % Define Grid
    xref = linspace( -1.0, 1.0, nx );
    yref = linspace( -1.0, 1.2, ny ); 

    % Interpolation Type (cubic or linear)
    interp_type = 'cubic';

    %% ndgrid Query Data 

    % Create Mesh with Ngrid or meshgrid
    [ X , Y ] = ndgrid(   xref, yref );

    % "Shifted" Query Points
    xq_shift = X+3.1*(xref(2)-xref(1));
    yq_shift = Y+3.1*(yref(2)-yref(1));

    % Create Test Surface
    Z = 3.0*peaks(X,Y);


    %% Meshgrid Query Data

    % meshgrid form
    [ Xm , Ym ] = meshgrid( xref, yref );

    % "Shifted" Query Points
    xq_shift_m = Xm+3.1*(xref(2)-xref(1));
    yq_shift_m = Ym+3.1*(yref(2)-yref(1));

    % Create Test Surface
    Zm = 3.0*peaks(Xm,Ym);


    %% mat-op-ex scheme

    % Get original size for reshaping
    nsize = size(xq_shift);

    % Get Sparse Matrix of interp2 operator
    [ M_shift ] = interp2_matrix( X, Y, xq_shift, yq_shift, interp_type );

    % Perform Test on Each Method
    t_mat_op_ex(1) = 0.0;
    for i=1:ni

        % New Z surface
        Z2 = Z + (i*0.3);

        % Shifted
        tstart = tic;
        Zval = reshape( M_shift*Z2(:), nsize );
        t_mat_op_ex(1) = t_mat_op_ex(1) + toc(tstart);


    end

    %% Interp2 - Must use meshgrid

    % Perform Test on Each Method
    t_interp2(1) = 0.0;
    for i=1:ni

        % New Z surface
        Z2 = Zm + (i*0.3);

        % Shifted
        tstart = tic;
        Zval2 = interp2( xref, yref, Z2, xq_shift_m, yq_shift_m, interp_type, NaN );
        t_interp2(1) = t_interp2(1) + toc(tstart);

    end


    %% Interpn 

    % Perform Test on Each Method
    t_interpn(1) = 0.0;
    for i=1:ni

        % New Z surface
        Z2 = Z + (i*0.3);

        % Shifted
        tstart = tic;
        Zval2 = interpn( xref, yref, Z2, xq_shift, yq_shift, interp_type, NaN );
        t_interpn(1) = t_interpn(1) + toc(tstart);

    end

    %% griddedInterpolant
        % Note: griddedInterpolant can store the X/Y mesh, but cannot
        % pre-calculate the query points, whereas interp2_matrix can. 

    % Create Function
    F = griddedInterpolant( X, Y, zeros(size(X)), interp_type, 'none' );

    % Perform Test on Each Method
    t_gridded(1) = 0.0;
    for i=1:ni

        % New Z surface
        Z2 = Z + (i*0.3);

        % Shifted
        tstart = tic;
        F.Values = Z2;
        Zval2 = F( xq_shift, yq_shift );
        t_gridded(1) = t_gridded(1) + toc(tstart);

        % Clear The F function, 
        F.Values(:,:) = 0.0;

    end

    %% GPU - interp2

    % Open device
    if (gpu_flag)

        % Reset GPU
        reset(gpu_device)

        % Create GPU Arrays
        Zg = gpuArray(Zm);
        xrefg = gpuArray(xref);
        yrefg = gpuArray(yref);
        xq_shift_g = gpuArray(xq_shift_m);
        yq_shift_g = gpuArray(yq_shift_m);

        % loop tests
        t_interp2_gpu(1)=0.0;
        for i=1:ni

            % New Z surface
            Z2 = Zg + (i*0.3);

            % Shifted
            tstart = tic;
            Zval2 = interp2( xrefg, yrefg, Z2, xq_shift_g, yq_shift_g, interp_type, NaN );
            t_interp2_gpu(1) = t_interp2_gpu(1) + toc(tstart);  

        end

        % Clear variables
        clearvars Zg xrefg yrefg xq_shift_g yq_shift_g  

    end

    %% GPU - mat-op-ex

    % Open device
    if (gpu_flag)

        % Create GPU Arrays
        M_shift_g = gpuArray(M_shift);
        Zg = gpuArray(Z);

        % loop tests
        t_mat_op_ex_gpu(1)=0.0;
        for i=1:ni

            % New Z surface
            Z2 = Zg + (i*0.3);

            % Shifted
            tstart = tic;
            Zval = reshape( M_shift_g*Z2(:), nsize );
            t_mat_op_ex_gpu(1) = t_mat_op_ex_gpu(1) + toc(tstart);

        end

        % Clear variables
        clearvars Zg M_shift_g 

    end
   
    % Store Times
    timer(1:4,itrial) =  [t_interpn ,t_interp2 ,t_gridded ,t_mat_op_ex ];
    if (gpu_flag)
        timer(5:6,itrial) = [t_interp2_gpu, t_mat_op_ex_gpu];
    else
        timer(5:6,itrial) = [nans(1), nans(1)];
    end
    
    % Normalize Times
    timer_relative(:,itrial) = timer(:,itrial)./ timer(1,itrial);
    
    % Print Status
    fprintf('\nTrial %i of %i',[itrial,ns]);
    
end

%% Plot Performance
% Open Figure
fig = figure(1); clf; hold on; grid on; 
fig.Units = 'inches';
fig.Position= [2,2,9,5];
   
% Title
title( { strcat( 'Time of Execution (Relative to interpn at each grid size) ',...
         sprintf(', %i Iterations',[ni] ) ),...
         ['CPU = ',cpu_name,', GPU = ',gpu_name]} )

% Plot cases
for i=1:6
    pp(i)=plot(nsa.*nsa,timer_relative(i,:),'.-',...
           'LineWidth',1.75,'MarkerSize',12);
end

% Set GPUs to dash
pp(5).LineStyle = '--';
pp(6).LineStyle = '--';

% Set Colors 
matop_color   = [0.1058 0.3411 0.502];
interp2_color = [0.478  0.294  0.294]; 
ginter_color  = [0.435  0.529  0.428];
interpn_color = [0.502  0.502  0.502];

% Recolor
pp(1).Color = interpn_color; 
pp(2).Color = interp2_color;
pp(3).Color = ginter_color; 
pp(4).Color = matop_color;
pp(5).Color = interp2_color; 
pp(6).Color = matop_color; 

% Get handle of the current axis
ax1 = gca;

% Set limits and scale of grid node axis
minx = nsa(1)*nsa(1);
maxx = nsa(end)*nsa(end);
xlim( [minx, maxx] )
set(ax1,'XScale','log')

% % Legend
legend(ax1,'interpn         - CPU',...
       'interp2         - CPU',...
       'grid.Interp.  - CPU',...
       'MAT-OP-EX  - CPU',...
       'interp2        - GPU',...
       'MAT-OP-EX  - GPU',...
       'Location','NorthWest')
   
% Label
xlabel(ax1,'Grid Nodes (Nx*Ny)')
ylabel(ax1,'Relative Time')

% Configure First axes for space for two more
pos = ax1.Position;
yshift = pos(1);
ax1.Position(2) = ax1.Position(2) + 2*yshift;
ax1.Position(4) = ax1.Position(4) - 2*yshift;

%% Create Second Axes - Grid Size in Mb
ax2 = axes('color','none');
ax2.Position = pos;
ax2.Position(2) = pos(2) + yshift;
ax2.Position(4) = 0.01;

% Compute estimates of data requirements for conventional interp2
if (strcmp(interp_type,'cubic')==1)
   mult = 16; % 4x4 stencil
else % linear, stencil is 4
   mult = 4;
end

% Compute Sizes (assuming matricies: Zq,Yq,F, linear X/Y)
minx_M = ( 3*8*minx*mult + 2*min_grid )/1e6;
maxx_M = ( 3*8*maxx*mult + 2*max_grid )/1e6;

% Define the 2nd X-axis
xlim([ minx_M, maxx_M])
xlabel(ax2,'Interp2/Interpn Memory Estimate (Mb)')
ax2.XColor = interp2_color;

%% Create Third Axes - Size of Mat-Op-Ex Operator
ax3 = axes('color','none');
ax3.Position = pos;
ax3.Position(4) = 0.01;

minx_M = minx_M + ( 16*minx*mult+8*minx )/1e6;
maxx_M = maxx_M + ( 16*maxx*mult+8*maxx )/1e6;

% Define the 2nd X-axis
xlim([ minx_M, maxx_M])
xlabel(ax3,'MAT-OP-EX Memory Estimate (Mb)')
ax3.XColor = matop_color;

 % Bold Everything
 set(findobj(gcf,'type','axes'),'FontName','Arial',...
'FontSize',10,'FontWeight','Bold');

