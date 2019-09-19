function M =  compute_neighborhood_casadi(position, r_comm, max_neig)
%% Build adjacency matrix
% Adjacency matrix M(i,j) (asymmetric in case of limited fov): 
% column j contains the neighbours i of agent j
% Note: agent j is not neig of itself (zeros on diagonal)
    
    import casadi.*

    N = length(position)/3;
    dist2_matrix = SX.sym('D',N,N);
    
    M = SX.zeros(N,N);
    
    for agent = 1:N
        for neig = 1:N
            if agent == neig
                dist2_matrix(agent,neig) = 1e10;
            elseif neig > agent
                dist2 = sum((position(3*(agent-1)+[1 2 3]) - ...
                    position(3*(neig-1)+[1 2 3])).^2);
                dist2_matrix(agent,neig) = dist2;
                dist2_matrix(neig,agent) = dist2; % for simmetry
%                 M(agent,neig) = if_else(dist2<r_comm^2, 1, M(agent,neig));
%                 M(neig,agent) = if_else(dist2<r_comm^2, 1, M(agent,neig));
            end
        end
    end
    % Vector counting the Number of Neighbours (nn)
    nn = sum(M,1);

    sorted_idx = SX.zeros(max_neig,N);
    for i = 1:N % for every column, order column
        % Check if constraint on topological distance is active
        %if nn(i) > max_neig
            % Sort the relative distance to neig
%             idx_dist2_to_sort = SX(1:N);
            to_sort = SX.ones(N,1);
            dist2_to_sort = dist2_matrix(:,i);
            for j = 1:max_neig % order only the first max_neig
                min_idx = 0;
                min_dist2 = SX(1e9); % high init value
                for k = 1:length(dist2_to_sort)
                    min_dist2 = if_else(dist2_to_sort(k) < min_dist2 & ...
                            dist2_to_sort(k) < r_comm^2 & ...
                            to_sort(k), ...
                            dist2_to_sort(k), min_dist2);
                    min_idx = if_else(dist2_to_sort(k) < min_dist2 & ...
                            dist2_to_sort(k) < r_comm^2 & ...
                            to_sort(k), ...
                            k, min_idx);
                   	to_sort(k) = if_else(dist2_to_sort(k) < min_dist2 & ...
                            dist2_to_sort(k) < r_comm^2 & ...
                            to_sort(k), ...
                            0, to_sort(k));
                    M(k,i) = if_else(dist2_to_sort(k) < min_dist2 & ...
                            dist2_to_sort(k) < r_comm^2 & ...
                            to_sort(k), ...
                            1, M(k,i));
                end
                
%                 idx_dist2_cell = horzsplit(idx_dist2_to_sort);
%                 default = 0;
%                 sorted_idx(j,i) = conditional(min_idx, idx_dist2_cell, default, false);
                sorted_idx(j,i) = min_idx;
%                 idx_dist2_to_sort = if_else(min_idx==0, [],...
%                     remove(idx_dist2_to_sort,min_idx-1,[]));
%                 dist2_to_sort = if_else(min_idx==0, [],...
%                     remove(dist2_to_sort,min_idx-1,[]));
            end

            % If it is active, cut the extra neighbors
%             to_keep = zeros(N,1);
%             cut_idx = find(sorted_idx(:,i),1,'last');
%             if cut_idx > 0 
%                 to_keep(sorted_idx(1:cut_idx,i)) = 1; 
%                 M(:,i) = to_keep;
%             end
        %end 
    end
    
end