% Forest map parameters
map.width = 50; % the forest is of size (width)x(width)
map.max_height = map.width; % maximum height of trees
map.nb_blocks = 1; % the number of blocks per row
map.street_width_perc = 0.8; % percentage of block that is empty

map.building_width = map.width/map.nb_blocks*(1-map.street_width_perc);
map.street_width = map.width/map.nb_blocks*map.street_width_perc;

map.building_shape = ('cylinder');

% Create buildings parameters
map = create_shifted_buildings(map);