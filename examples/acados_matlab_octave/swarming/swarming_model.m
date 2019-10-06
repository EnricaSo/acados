function model = swarming_model(S)

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
u_ref = S.u_migration; % reference direction of velocity for all agents
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

%% Nonlinear least squares

% Weights
W_sep = 1; 
W_dir = 1;
W_nav = 2;
W_u = 2e-1; % Penalization of high values of the control input variables

sym_sep = SX.zeros(N*(N-1),1);
sym_dir = SX.zeros(N,1);
sym_nav = SX.zeros(N,1);

%ny = N*(N+1);

% Neighborhood matrix
M = ones(N,N) - eye(N,N);
% M = compute_neighborhood_casadi(pos, r_comm, max_neig);

% For every agent define the nonlinear_ls terms
for agent = 1:N
    
    % Get the index triplet related to the current agent
    agent_idx = [1,2,3]' + 3*(agent-1)*ones(3,1);
    
    % For every neighbor, compute the distance to the current agent
    for neig = 1:(N-1)
        if neig < agent
            neig_idx = [1,2,3]' + 3*(neig-1)*ones(3,1);
        else
            neig_idx = [1,2,3]' + 3*(neig)*ones(3,1);
        end
        % Separation term
        pos_rel = pos(neig_idx)-pos(agent_idx);
        sym_sep((agent-1)*(N-1)+neig) = 1/(N-1)*(pos_rel'*pos_rel - d_ref^2);
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

%% Constraints

sym_dist = SX.zeros(N*(N-1)/2,1);

k = 1;
for agent = 1:(N-1)
    agent_idx = [1,2,3]' + 3*(agent-1)*ones(3,1);
    for neig = (agent+1):N
        neig_idx = [1,2,3]' + 3*(neig-1)*ones(3,1);
        pos_rel = pos(neig_idx)-pos(agent_idx);
        sym_dist(k) = pos_rel(1)^2+pos_rel(2)^2+pos_rel(3)^2;
        k = k+1;
    end
end

% Constraints on control inputs and distances
expr_h = vertcat(sym_u, ...
                 sym_dist); 
% expr_h_e = sym_x;

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

