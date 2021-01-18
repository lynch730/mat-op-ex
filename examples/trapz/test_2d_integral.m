
close all
clear; clc;
addpath('../../source/')

% This script demonstrates the matrix extraction method for solving a
% spherical integration problem with symmetry along the phi (azimuth) axis.
% (effectively a 2D problem)

%% Reference Data Setup

% Get reference data
nt = 18; % theta (X)
nr = 12; % radius (Y) 

% Define Grid Limits
theta(:,1)  = linspace( 0.0, pi, nt ); % X
radius(:,1) = linspace( 0.0, 5.0, nr ); % Y

% Create Mesh with Ngrid or meshgrid
[ T, R ] = ndgrid( theta, radius);

% differential elements
V_T = 2.0.*pi.*sin(theta); % phi-integrated differential surface
V_R = radius.*radius;      % Radial factor

%% Matrix Extraction

% Simpson Integrate over Theta
[ M_theta_simp ] = trapz_matrix2( T, R,  1, 'simpsons' , V_T);

% Simpson Integrate over R
[ M_radius_simp ] = trapz_matrix1( radius, 'simpsons', V_R );

% Trapz Integrate over Theta
[ M_theta_trapz ] = trapz_matrix2( T, R,  1, 'trapz' , V_T);

% TrapZ Integrate over R
[ M_radius_trapz ] = trapz_matrix1( radius, 'trapz', V_R );

%% Application 

% Create Test Surface of R
Z = R(:) ; 

% Linearize
Zlin = Z(:);

% Solve
integral_simps = M_radius_simp*M_theta_simp*Zlin;
integral_trapz = M_radius_trapz*M_theta_trapz*Zlin;

% Error 
val = 625*pi; % known value, for R=5;
error_simps = 100.0*abs( (integral_simps - val) ./ val); 
error_trapz = 100.0*abs( (integral_trapz - val) ./ val); 

% Display
fprintf('Spherical Integration %% Deviation');
fprintf('\n Trapz    | %12.10f   ',error_trapz);
fprintf('\n Simpsons | %12.10f \n',error_simps);

