function check = check_reformulation(model, gnsf, print_info)
%
%   This file is part of acados.
%
%   acados is free software; you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation; either
%   version 3 of the License, or (at your option) any later version.
%
%   acados is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with acados; if not, write to the Free Software Foundation,
%   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%   Author: Jonathan Frey: jonathanpaulfrey(at)gmail.com

%% Description:
% this function takes the implicit ODE/ index-1 DAE and a gnsf structure
% to evaluate both models at num_eval random points x0, x0dot, z0, u0;
% if for all points the relative error is <= TOL, the function will return
% 1, otherwise it will give an error.

import casadi.*

TOL = 1e-14;
num_eval = 10;

% get dimensions
nx  = gnsf.nx;
nu  = gnsf.nu;
nz  = gnsf.nz;

nx1 = gnsf.nx1;
nx2 = gnsf.nx2;
% n_out = gnsf.n_out;
% ny = gnsf.ny;
% nuhat = gnsf.nuhat;


% get model matrices
A  = gnsf.A;
B  = gnsf.B;
C  = gnsf.C;
E  = gnsf.E;
c  = gnsf.c;

L_x    = gnsf.L_x;
L_xdot = gnsf.L_xdot;
L_z    = gnsf.L_z;
L_u    = gnsf.L_u;

A_LO = gnsf.A_LO;

I_x1 = 1:nx1;
I_x2 = nx1+1:nx;


% get casadi variables
x = gnsf.x;
xdot = gnsf.xdot;
z = gnsf.z;
u = gnsf.u;
y = gnsf.y;
uhat = gnsf.uhat;

% create functions
impl_ode_fun = Function(['impl_ode_fun'], {x, xdot, u, z}, {model.f_impl_expr});
phi_fun = Function('phi_fun',{y,uhat}, {gnsf.phi_expr});
f_lo_fun = Function('f_lo_fun',{x(1:nx1), xdot(1:nx1), z, u}, {gnsf.f_lo_expr});



for i_check = 1:num_eval
    
    % generate random values
    x0    = rand(nx, 1);
    x0dot = rand(nx, 1);
    z0    = rand(nz, 1);
    u0    = rand(nu, 1);
    
    
    % eval f_impl;
    f_impl_val = full(impl_ode_fun(x0, x0dot, u0, z0));

    % eval gnsf
    y0 = L_x * x0(I_x1) + L_xdot * x0dot(I_x1) + L_z * z0;
    uhat0 = L_u * u0;

    gnsf_val1 = (A * x0(I_x1) + B * u0 + ...
        C * phi_fun( y0, uhat0) + c) - E * [x0dot(I_x1); z0];
    
    if nx2 > 0 % eval LOS
%         gnsf_val2 = A_LO * x0(I_x2) + ...
%             gnsf.f_lo_fun(x0(I_x1), x0dot(I_x1), z0, u0) - x0dot(I_x2);
        gnsf_val2 =  A_LO * x0(I_x2) + ...
            f_lo_fun(x0(I_x1), x0dot(I_x1), z0, u0) - x0dot(I_x2);
        gnsf_val = full([gnsf_val1; gnsf_val2 ]);
    else
        gnsf_val = full(gnsf_val1);
    end
    
    % compute error and check
    rel_error = norm(f_impl_val - gnsf_val) / norm(f_impl_val);

    if rel_error > TOL
        abs_error = gnsf_val - f_impl_val;
        T = table(f_impl_val, gnsf_val, abs_error);
        disp(T)
        return
        error('transcription failed; rel_error > TOL');
        check = 0;
    end            
end

if print_info
    disp(['model reformulation checked: relative error <= TOL = ', num2str(TOL)]);
end
check = 1;
    
%% helpful for debugging
% % use in calling function and compare
% % compare f_impl(i) with gnsf_val1(i);
% 

%     nx  = gnsf.nx;
%     nu  = gnsf.nu;
%     nz  = gnsf.nz;
%     nx1 = gnsf.nx1;
%     nx2 = gnsf.nx2;
%     
%         A  = gnsf.A;
%     B  = gnsf.B;
%     C  = gnsf.C;
%     E  = gnsf.E;
%     c  = gnsf.c;
% 
%     L_x    = gnsf.L_x;
%     L_z    = gnsf.L_z;
%     L_xdot = gnsf.L_xdot;
%     L_u    = gnsf.L_u;
%     
%     A_LO = gnsf.A_LO;
% 
%     x0 = rand(nx, 1);
%     x0dot = rand(nx, 1);
%     z0 = rand(nz, 1);
%     u0 = rand(nu, 1);
%     I_x1 = 1:nx1;
%     I_x2 = nx1+1:nx;
%         
%     y0 = L_x * x0(I_x1) + L_xdot * x0dot(I_x1) + L_z * z0;
%     uhat0 = L_u * u0;
% 
%     gnsf_val1 = (A * x(I_x1) + B * u + ...
%         C * phi_current + c) - E * [xdot(I_x1); z];
%     gnsf_val1 = gnsf_val1.simplify();
%     
% %     gnsf_val2 = A_LO * x(I_x2) + gnsf.f_lo_fun(x(I_x1), xdot(I_x1), z, u) - xdot(I_x2);
%     gnsf_val2 =  A_LO * x(I_x2) + gnsf.f_lo_fun(x(I_x1), xdot(I_x1), z, u) - xdot(I_x2);
% 
%         
%     gnsf_val = [gnsf_val1; gnsf_val2];
%     gnsf_val = gnsf_val.simplify();
%     f_impl_expr = f_impl_expr.simplify();
% keyboard


end