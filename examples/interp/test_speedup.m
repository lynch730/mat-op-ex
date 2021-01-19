
close all
clear; clc;
addpath('../../source/')

% This script tests the speedup of interp2_matrix for 2D matrix
% interpolation on a fixed grid, relative to interp2, griddedInterpolant,
% and with a GPU, if available. Two sets of queries are tested, one where
% query points are shifted by uniform amounts from the reference grid ("shift"), and
% one where they are more randomized ("scatter")

%% Hardcode CPU/GPU info
cpu_name = 'Intel i7-5820K';
gpu_name = 'NVIDIA GTX 1070';

% Number of iterations to time
ni = 1000;

% Loop number of grid points trials
ns = 20;
nsa = round(linspace(10,600,ns));

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
title( { 'Time of Execution (Relative to interpn at each grid size) ',...
         sprintf('%i Iterations',[ni] ),...
         ['CPU = ',cpu_name,', GPU = ',gpu_name]} )

% Plot CPU's
for i=1:4
    pp(i)=plot(nsa.*nsa,timer_relative(i,:),'.-',...
           'LineWidth',1.5,'MarkerSize',15);
end

% Plot GPU's
for i=5:6
    pp(i)=plot(nsa.*nsa,timer_relative(i,:),'.--',...
           'LineWidth',1.5,'MarkerSize',15);
end

xlim( [min(nsa.*nsa), max(nsa.*nsa)] )

% % Legend
legend('interpn         - CPU',...
       'interp2         - CPU',...
       'grid.Interp.  - CPU',...
       'MAT-OP-EX  - CPU',...
       'interp2        - GPU',...
       'MAT-OP-EX  - GPU')
   
% Label
xlabel('Grid Nodes (Nx*Ny)')
ylabel('Relative Time')
set(gca,'XScale','log')
     
 set(findobj(gcf,'type','axes'),'FontName','Arial',...
'FontSize',10,'FontWeight','Bold');



