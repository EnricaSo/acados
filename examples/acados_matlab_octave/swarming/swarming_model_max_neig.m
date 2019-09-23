function model = swarming_model_neig(S)

% SWARMING_MODEL - Function that describes the dynamics of the swarm and
% the cost function for controlling it.
%
% Swarming  is the behavior of collective and ordered motion. It can be 
% obtained through the combination of the following rules:
%
% Separation rule: drives the agents to a reference inter-agent ...
%       distance (d_ref)
% Direction rule: make the agents' velocity converge to a ...
%       reference direction (u_ref)
% Navigation rule: make the agents' speed converge to a reference ...
%       value (v_ref)
%

import casadi.*

%% Add parameters

%% Rename swarming parameters

N = S.nb_agents; % nb of agents
max_neig = S.max_neig; % number of neighbours
v_ref = S.v_flock; % reference speed for all agents
u_ref = S.u_migration; % reference velocity direction for all agents (unit vector)
d_ref = S.d;  % reference distance among every couple of neighboring agents
max_a = S.max_a;
r_comm = S.r;

%% System dimensions

nx = 6 * N; % nb of state variables
nu = 3 * N; % nb of control inputs

%% Named symbolic variables

pos = SX.sym('p', 3*N); % 3D positions of the agents [m]
vel = SX.sym('v', 3*N); % 3D velocities of the agents [m/s]
u = SX.sym('a', 3*N); % 3D acceleration to apply to agents [m/s^2]

%% Unnamed symbolic variables

sym_x = vertcat(pos, vel);
sym_xdot = SX.sym('xdot', nx, 1);
sym_u = u;

%% Dynamics

expr_f_expl = vertcat(vel, ...
                      u);
expr_f_impl = expr_f_expl - sym_xdot;

%% Constraints

expr_h = sym_u; % constraints only on control inputs, for now
% expr_h_e = sym_x;

%% Nonlinear least squares

% Weights
W_sep = 1; 
W_dir = 1;
W_nav = 2;
W_u = 1e-1; % Penalization of high values of the control input variables

sym_sep = SX.zeros(N*(N-1),1);
sym_dir = SX.zeros(N,1);
sym_nav = SX.zeros(N,1);

%ny = N*(N+1);

% Neighborhood matrix
% M = ones(N,N) - eye(N,N);
[~, sorted_neig] = compute_neighborhood_casadi(pos, r_comm, max_neig);

% For every agent define the nonlinear_ls terms
for agent = 1:N
    
    % Get the index triplet related to the current agent
    agent_idx = [1,2,3]' + 3*(agent-1)*ones(3,1);
    neigs = sorted_neig(:,agent);
    % For every neighbor, compute the distance to the current agent
    for j = 1:max_neig
        neig = neigs(j);
        neig_idx = [1,2,3]' + 3*(neig)*ones(3,1);
        % Separation term
        pos_rel_cell = vertsplit(pos-repmat(pos(agent_idx),N,1));
        pos_rel_default = [d_ref;0;0];
        pos_rel = SX.zeros(3,1);
        neig_idx_x = neig_idx(1)-1;
        neig_idx_y = neig_idx(2)-1;
        neig_idx_z = neig_idx(3)-1;
        pos_rel(1) = conditional(neig_idx_x, pos_rel_cell, pos_rel_default(1), false);
        pos_rel(2) = conditional(neig_idx_y, pos_rel_cell, pos_rel_default(2), false);
        pos_rel(3) = conditional(neig_idx_z, pos_rel_cell, pos_rel_default(3), false);
        sym_sep((agent-1)*(N-1)+j) = 1/(N-1)*(pos_rel'*pos_rel - d_ref^2);
    end
    vel_agent = vel(agent_idx);
    % Direction term
    sym_dir(agent) = 1 - (vel_agent'*u_ref)^2/(vel_agent'*vel_agent);
    % Navigation term
    sym_nav(agent) = vel_agent'*vel_agent - v_ref^2;
end

sym_sep = W_sep * sym_sep;
sym_dir = W_dir * sym_dir;
sym_nav = W_nav * sym_nav;

% Assemble expr_y
expr_y = vertcat(sym_sep, sym_dir, sym_nav, W_u*sym_u);
expr_y_e = vertcat(sym_sep, sym_dir, sym_nav);

ny = length(expr_y);
ny_e = length(expr_y_e);

%% Populate structure

model.nx = nx;
model.nu = nu;
model.ny = ny;
model.ny_e = ny_e;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;
model.expr_h = expr_h;
% model.expr_h_e = expr_h_e;
model.expr_y = expr_y;
model.expr_y_e = expr_y_e;

