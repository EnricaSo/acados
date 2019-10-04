%% Arguments

% Import swarming param, in S structure
run('/home/esoria/Developer/sp2018_uavSim/swarming/param/param_swarm');

% Overwrite and add some param
% S.nb_agents = 5;
% S.max_neig = 2;
S.max_a = 2;
S.d = 5;

% Rename param
N = S.nb_agents; % nb of agents
max_neig = S.max_neig; % number of neighbours
v_ref = S.v_flock;
u_ref = S.u_migration;
d_ref = S.d;
max_a = S.max_a;
r_coll = S.r_coll;

%% Call function
model = swarming_model_max_neig(S);