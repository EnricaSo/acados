function cost = compute_cost_offline(S, model, pos_history, vel_history, u_history)
% compute_cost_offline - function that allows to retrieve a vector of
% values with the MPC cost at every time step.

%% Rename parameters

% Swarming
N = S.nb_agents; % nb of agents
max_neig = S.max_neig; % number of neighbours
v_ref = S.v_flock; % reference speed for all agents
u_ref = S.u_migration; % reference direction of velocity for all agents
d_ref = S.d;  % reference distance among every couple of neighboring agents
max_a = S.max_a;
r_comm = S.r;
r_coll = S.r_coll;

% Weights
W_sep = model.W_sep;
W_dir = model.W_dir;
W_nav = model.W_nav;
W_u = model.W_u; % Penalization of high values of the control input variables


%% Main loop to compute the MPC cost

[nb_steps,~] = size(pos_history(1:(end-1),:));
cost = zeros(nb_steps,1);

for step = 1:nb_steps
    
    % Neighborhood matrix
    M = ones(N,N) - eye(N,N);
    % M = compute_closest_neig(pos, r_comm, max_neig);
    
    sym_sep = zeros(N*(N-1),1);
    sym_dir = zeros(N,1);
    sym_nav = zeros(N,1);
    
    pos = (pos_history(step,:))';
    vel = (vel_history(step,:))';
    u = (u_history(step,:))';
    
    % For every agent define the nonlinear_ls terms
    for agent = 1:N
        
        % Get the index triplet related to the current agent
        agent_idx = [1,2,3]' + 3*(agent-1)*ones(3,1);
        
        % For every neighbor, compute the distance to the current agent
        for j = 1:(N-1)
            if j < agent
                neig = j;
            else
                neig = j+1;
            end
            neig_idx = [1,2,3]' + 3*(neig-1)*ones(3,1);
            % Separation term
            pos_rel = pos(neig_idx)-pos(agent_idx);
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
    cost(step) = sum(sym_sep.^2) + sum(sym_dir.^2) + sum(sym_nav.^2) + W_u*sum(u.^2);
    
end

end