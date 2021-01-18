
close all
clear; clc;
addpath('core/');

% This script demonstrates the matrix extraction method of solving trapz
% using matrix methods, applied to a problem of 1D grids. 

%% Reference Data Setup

% Get reference data
Nx = 15;

% Define Grid Limits
xref(:,1) = linspace( -2.0, 3.0, Nx );

% Create Test Surface
Z(:,1) = 10.0.*xref.*xref + 0.8.*xref + 1.2;

%% Test Integration
   % Test to show the trapz and Simpsons methods work

% Trapz Test with Matrix Method
Int_trap = trapz_matrix1( xref , 'trapz' )*Z;

% Composite Simpsons Rule with Matrix Method
Int_simp = trapz_matrix1( xref , 'simpsons' )*Z;

% Exact solution
Int_test(:,1) = 124.0+2.0/3.0;

% Deviations
error_trap = 100.0.*abs( (Int_trap-Int_test)./Int_test );
error_simp = 100.0.*abs( (Int_simp-Int_test)./Int_test );

%% Test V vector
   % Test to show the V mult Array Works

% Test Scalar Multiplier
Int_vmult1 = trapz_matrix1( xref , 'simpsons' , 5.0 )*Z ./ 5.0;

% Test Vector Multiplier
Int_vmult2 = trapz_matrix1( xref , 'simpsons' , ones(Nx,1)*5.0 )*Z ./5.0;

% Deviations
error_vscalar = 100.0.*abs( (Int_vmult1-Int_simp)./Int_simp );
error_varray  = 100.0.*abs( (Int_vmult2-Int_simp)./Int_simp );

%% Display Deviations
fprintf('trapz_matrix1, %% Deviation');
fprintf('\n Trapz    | %16.10e   ',error_trap);
fprintf('\n Simpsons | %16.10e   ',error_simp);
fprintf('\n Scalar V | %16.10e   ',error_vscalar);
fprintf('\n Vector V | %16.10e \n',error_varray);

