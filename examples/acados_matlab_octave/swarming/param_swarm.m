% Variables to be set
S.is_active_migration = false;
S.is_active_goal = true;
S.is_active_arena = false;
S.is_active_spheres = false;
S.is_active_cyl = true;
S.draw_plot = false;
S.draw_state_var = false;
S.draw_neig_links = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of agents
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S.nb_agents = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Max radius of influence - Metric distance
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S.r = 150;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Max number of neighbors - Topological distance
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S.max_neig = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Max field of view
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P.aov = 120;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radius of collision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S.r_coll = 0.5;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cylindric obstacles parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('map','var') & ~isempty(map))

    nb_obstacles = length(map.buildings_east);
    cylinder_radius = map.building_width / 2;

    S.cylinders = [
        map.buildings_north'; % x_obstacle
        map.buildings_east'; % y_obstacle
        repmat(cylinder_radius, 1, nb_obstacles)]; % r_obstacle
    S.n_cyl = length(S.cylinders(1, :));
else
    S.cylinders = 0;
    S.n_cyl = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Speed parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Migration direction
S.u_migration = [1 0 0]';

% Migration speed
S.v_flock = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity and acceleration bounds for the agents
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S.max_a = sqrt(3*4);
S.max_v = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial position and velocity for the swarm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S.P0 	= [-10,20,20]'; %                                [m]
S.P     = 10; %                                          [m]
S.V0 	= [0,0,0]'; %                                    [m/s]
S.V 	= 0; %                                           [m/s]
S.seed = 5;

rng(S.seed);
S.Pos0      = S.P0 + S.P * rand(3,S.nb_agents);
S.Vel0      = S.V0 + S.V * rand(3,S.nb_agents);
