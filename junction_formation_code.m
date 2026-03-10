close all
clearvars

tic

%This code randomly chooses connection point amongst all points within the
%given threshold, weighted by inverse distance (more likely to connect at
%the closest points). Diffusion of rods from Broersma theory is used for
%translational and rotational diffusion. Measures the interior angles of 
%L-junctions from 0-180 degrees (X- and T-junctions are measured obtusely, 
%too, and forced to be acute in post processing)

% Code written by Brittany Bannish, Key Crosley, and Marcus Case. Last
% updated 3/10/2026


%% Parameters
L = 115; % side length of the cube in um
segment_length = 6; %12; % length of each line segment in um
num_segments = 200; % number of line segments
Tend = 420; %length of time in seconds to run simulation for
dt = 0.01; % time step in s
num_steps = ceil(Tend/dt); % number of steps in the simulation
%rotational_diffusion = 0.2; % scale factor of the random orientation (rotational) movement
%translational_diffusion = 1; %scale factor of the random translational movement
proximity_threshold = 0.1; % threshold distance for segments to stick together in um
discrete_size = 20; % number of discretization points per segment
plot_interval = 1000;

%Parameters for diffusion coefficients
kB = 1.380649*10^(-23); %Boltzmann's constant (m^2 * kg) / (s^2 * K)
T = 310.15; %human body temperature in K (room temp, which was experiments, is about 296)
eta = 1.8*6.193*10^(-4); %dynamic viscosity of blood plasma is 1.8-times the viscosity of water in kg / (m * s)
drod = 1*10^(-7); %diameter of rod in m (=100 nm)
Lrod = segment_length*10^(-6); %length of rod in m

% Calculate log term (used in all formulas)
log_term = log(Lrod/drod);
log_2term = log(2*Lrod/drod);

% Broersma translational diffusion coefficient - PARALLEL to rod axis
D_parallel = kB*T*(log_term - 0.114 - 0.15/(log_2term) - 13.5/(log_2term^2) + 37/(log_2term^3) - 22/(log_2term^4))/(2*pi*eta*Lrod);
D_parallel = D_parallel * 10^12; % Convert to um^2/s

% Broersma translational diffusion coefficient - PERPENDICULAR to rod axis
D_perp = kB*T*(log_term + 0.866 - 0.15/log_2term - 8.1/(log_2term^2) + 18/(log_2term^3) - 9/(log_2term^4))/(4*pi*eta*Lrod);
D_perp = D_perp * 10^12; % Convert to um^2/s

% Calculate rotational diffusion coefficient - Broersma corrected (you already have this)
D_rot = 3*kB*T*(log_term - 0.446 - 0.2/log_2term - 16/(log_2term^2) + 63/(log_2term^3) - 62/(log_2term^4))/(pi*eta*Lrod^3); %measured in rad^2/s

%parameters for nearness calculation
nearness = [];
near_proximity = 0.25;
countnear = 0;

%rng(1) % Gives the same random numbers every time, comment out if needed

% Video
videoFileName = 'run1.mp4';
v = VideoWriter(videoFileName, 'MPEG-4');
v.FrameRate = 5; % Adjust frame rate as needed
open(v);
countstep=0;

%% Initialize the positions of the endpoints of the line segments
% Generate random centers for segments with margin for segment length
margin = segment_length / 2;
centers = rand(num_segments, 3) * (L - 2*margin) + margin;

% Generate random directions
directions = randn(num_segments, 3); % random initial directions
directions = directions ./ vecnorm(directions, 2, 2); % normalize directions

% Calculate endpoints from centers
p1 = centers - (segment_length/2) * directions;
p2 = centers + (segment_length/2) * directions;

% Initialize cell array for discretized segments
P = cell(1, num_segments);
% Discretize each segment
for i = 1:num_segments
    P{i} = [linspace(p1(i,1), p2(i,1), discrete_size);
        linspace(p1(i,2), p2(i,2), discrete_size);
        linspace(p1(i,3), p2(i,3), discrete_size)];
end

% Initialize connection tracking - stores [seg1, seg2, point1, point2] for each connection
connections = [];
junction_type = [];
theta = [];
junction_formation_time = [];  % Track when each junction forms
junction_points = [];  % Track which discretized points are involved [seg1_point, seg2_point]

Lform = [];
count_L = 0; % counts number of L junctions
seg_in_comp = zeros(num_segments,1); % initialize array of segments that are involved in a connection


%% Plotting setup
figure(1);
view(3);
hold on;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
xlim([0 L]); ylim([0 L]); zlim([0 L]);
title('Fiber Diffusion and Junction Formation in 3D');
colors = lines(num_segments);

% Initialize plot handles
segment_plots = gobjects(num_segments, 1);
for i = 1:num_segments
    segment_plots(i) = plot3(P{i}(1,:), P{i}(2,:), P{i}(3,:), '-', ...
        'Color', colors(i,:), 'LineWidth', 2);
end
drawnow;


connected = false(num_segments, num_segments);


%% Main simulation loop (num_steps = 42000 for a 420 w 0.01 timesteps
for step = 1:num_steps

    % Build connected_segments from connections matrix for component identification
    connected_segments = cell(1, num_segments);
    for i = 1:num_segments
        connected_segments{i} = [];
    end

    % Populate connected_segments from connections matrix
    for conn_idx = 1:size(connections, 1)
        seg1 = connections(conn_idx, 1);
        seg2 = connections(conn_idx, 2);
        connected_segments{seg1} = [connected_segments{seg1}, seg2];
        connected_segments{seg2} = [connected_segments{seg2}, seg1];
    end

    % Identify connected components (groups of connected segments)
    visited = false(1, num_segments);
    components = {};
    for i = 1:num_segments
        if ~visited(i)
            % Start new component with segment i
            component = [];
            stack = i;

            while ~isempty(stack)
                current = stack(end);
                stack(end) = [];

                if ~visited(current)
                    visited(current) = true;
                    component = [component, current];

                    % Add all connected segments to stack
                    stack = [stack, connected_segments{current}];
                end
            end

            components{end+1} = component;
        end
    end

    % Move each component as a rigid body
    for comp_idx = 1:length(components)
        component = components{comp_idx};

        % Generate random movement and rotation for this component.
        % translational noise will be calculated separately for the single
        % segments and connected components in the loop below
        noise_rot = randn(1,3)*sqrt(2*D_rot*dt);


        if length(component) == 1
            % Single segment - move with anisotropic diffusion
            i = component(1);
            center = mean(P{i}, 2);
            current_direction = P{i}(:,end) - P{i}(:,1);
            current_direction = current_direction / norm(current_direction);

            % Generate random displacement components
            noise_parallel = randn(1,1)*sqrt(2*D_parallel*dt);  % Along rod axis
            noise_perp = randn(2,1)*sqrt(2*D_perp*dt);          % Perpendicular to rod axis

            % Create perpendicular basis vectors to rod axis
            % First perpendicular vector
            if abs(current_direction(1)) < 0.9
                perp1 = cross(current_direction, [1;0;0]);
            else
                perp1 = cross(current_direction, [0;1;0]);
            end
            perp1 = perp1 / norm(perp1);

            % Second perpendicular vector
            perp2 = cross(current_direction, perp1);
            perp2 = perp2 / norm(perp2);

            % Combine displacements in rod reference frame
            displacement = noise_parallel * current_direction + ...
                noise_perp(1) * perp1 + ...
                noise_perp(2) * perp2;

            % Apply random movement to center
            new_center = center + displacement;

            % Apply random orientation change (this part stays the same)
            axis = noise_rot';
            angle = norm(axis);
            if angle > 0
                axis = axis / angle;
                new_direction = current_direction + angle * cross(axis, current_direction);
                new_direction = new_direction / norm(new_direction);
            else
                new_direction = current_direction;
            end

            % Update segment endpoints
            p1_new = new_center - (segment_length/2) * new_direction;
            p2_new = new_center + (segment_length/2) * new_direction;

            % Apply boundary conditions
            p1_new = max(0, min(L, p1_new));
            p2_new = max(0, min(L, p2_new));

            % Ensure segment length is maintained
            actual_direction = p2_new - p1_new;
            actual_length = norm(actual_direction);
            if actual_length > 0
                actual_direction = actual_direction / actual_length;
                p2_new = p1_new + segment_length * actual_direction;
            end

            % =============== Re-discretize the segment
            P{i} = [linspace(p1_new(1), p2_new(1), discrete_size);
                linspace(p1_new(2), p2_new(2), discrete_size);
                linspace(p1_new(3), p2_new(3), discrete_size)];

        else
            % Connected component - move as rigid body while preserving connections.
            % use averaged diffusion coefficients

            D_avg = (D_parallel + 2*D_perp)/3;  % Average of parallel and perpendicular
            noise_trans = randn(1,3)*sqrt(2*D_avg*dt);

            % Store original connection points before any movement
            original_connection_points = [];
            for conn_idx = 1:size(connections, 1)
                seg1 = connections(conn_idx, 1);
                seg2 = connections(conn_idx, 2);
                point1_idx = connections(conn_idx, 3);
                point2_idx = connections(conn_idx, 4);

                % Only store connections within this component
                if ismember(seg1, component) && ismember(seg2, component)
                    original_connection_points(end+1,:) = [conn_idx, seg1, seg2, point1_idx, point2_idx];
                end
            end

            % Calculate center of mass of the component
            all_points = [];
            for i = component
                all_points = [all_points, P{i}];
            end
            com = mean(all_points, 2);

            % % Apply translation
            % new_com = com + (1/length(components{comp_idx}))*noise_trans'; %have bigger components move slower
            % for i = component
            %     P{i} = P{i} + (1/length(components{comp_idx}))*noise_trans'; %have bigger components move slower
            % end
            new_com = com + (1/(length(components{comp_idx}))^3)*noise_trans'; %have bigger components move MUCH slower
            for i = component
                P{i} = P{i} + (1/(length(components{comp_idx}))^3)*noise_trans'; %have bigger components move MUCH slower
            end

            % Apply rotation around center of mass
            %axis = (1/length(components{comp_idx}))*noise_rot'; %have bigger components move slower
            axis = (1/(length(components{comp_idx}))^3)*noise_rot'; %have bigger components move MUCH slower
            angle = norm(axis);
            if angle > 0
                axis = axis / angle;

                % Create rotation matrix for small angles
                K = [0, -axis(3), axis(2);
                    axis(3), 0, -axis(1);
                    -axis(2), axis(1), 0];
                R = eye(3) + sin(angle) * K + (1-cos(angle)) * K^2;

                % Apply rotation to each segment around component
                for i = component
                    P{i} = R * (P{i} - new_com) + new_com;
                end
            end

            % Enforce connection constraints after movement
            for conn_row = 1:size(original_connection_points, 1)
                conn_idx = original_connection_points(conn_row, 1);
                seg1 = original_connection_points(conn_row, 2);
                seg2 = original_connection_points(conn_row, 3);
                point1_idx = original_connection_points(conn_row, 4);
                point2_idx = original_connection_points(conn_row, 5);

                % Get current positions of connection points
                point1 = P{seg1}(:, point1_idx);
                point2 = P{seg2}(:, point2_idx);

                % Calculate midpoint
                midpoint = (point1 + point2) / 2;

                % Move both connection points to the midpoint
                displacement1 = midpoint - point1;
                displacement2 = midpoint - point2;

                % Apply constraint forces to maintain segment lengths
                % For segment 1
                p1_seg1 = P{seg1}(:, 1);
                p2_seg1 = P{seg1}(:, end);
                center1 = (p1_seg1 + p2_seg1) / 2;

                % Move segment 1 to satisfy connection constraint
                P{seg1} = P{seg1} + displacement1;

                % Ensure segment 1 maintains its length
                new_p1_seg1 = P{seg1}(:, 1);
                new_p2_seg1 = P{seg1}(:, end);
                current_dir1 = new_p2_seg1 - new_p1_seg1;
                current_length1 = norm(current_dir1);
                if current_length1 > 0
                    current_dir1 = current_dir1 / current_length1;
                    new_center1 = (new_p1_seg1 + new_p2_seg1) / 2;
                    new_p1_seg1 = new_center1 - (segment_length/2) * current_dir1;
                    new_p2_seg1 = new_center1 + (segment_length/2) * current_dir1;

                    % =============== Re-discretize segment 1
                    P{seg1} = [linspace(new_p1_seg1(1), new_p2_seg1(1), discrete_size);
                        linspace(new_p1_seg1(2), new_p2_seg1(2), discrete_size);
                        linspace(new_p1_seg1(3), new_p2_seg1(3), discrete_size)];
                end

                % For segment 2
                p1_seg2 = P{seg2}(:, 1);
                p2_seg2 = P{seg2}(:, end);
                center2 = (p1_seg2 + p2_seg2) / 2;

                % Move segment 2 to satisfy connection constraint
                P{seg2} = P{seg2} + displacement2;

                % Ensure segment 2 maintains its length
                new_p1_seg2 = P{seg2}(:, 1);
                new_p2_seg2 = P{seg2}(:, end);
                current_dir2 = new_p2_seg2 - new_p1_seg2;
                current_length2 = norm(current_dir2);
                if current_length2 > 0
                    current_dir2 = current_dir2 / current_length2;
                    new_center2 = (new_p1_seg2 + new_p2_seg2) / 2;
                    new_p1_seg2 = new_center2 - (segment_length/2) * current_dir2;
                    new_p2_seg2 = new_center2 + (segment_length/2) * current_dir2;

                    % =========== Re-discretize segment 2
                    P{seg2} = [linspace(new_p1_seg2(1), new_p2_seg2(1), discrete_size);
                        linspace(new_p1_seg2(2), new_p2_seg2(2), discrete_size);
                        linspace(new_p1_seg2(3), new_p2_seg2(3), discrete_size)];
                end
            end

            % Apply boundary conditions to entire component
            for i = component
                % Get current endpoints
                p1_current = P{i}(:,1);
                p2_current = P{i}(:,end);

                % Apply boundary reflection to endpoints
                p1_new = p1_current;
                p2_new = p2_current;

                % Reflect if outside boundaries
                for dim = 1:3
                    if p1_current(dim) < 0
                        p1_new(dim) = -p1_current(dim);
                    elseif p1_current(dim) > L
                        p1_new(dim) = 2*L - p1_current(dim);
                    end

                    if p2_current(dim) < 0
                        p2_new(dim) = -p2_current(dim);
                    elseif p2_current(dim) > L
                        p2_new(dim) = 2*L - p2_current(dim);
                    end
                end

                % Maintain segment length
                direction = p2_new - p1_new;
                current_length = norm(direction);
                if current_length > 0
                    direction = direction / current_length;
                    p2_new = p1_new + segment_length * direction;
                end

                % =========== Re-discretize to maintain straight line
                P{i} = [linspace(p1_new(1), p2_new(1), discrete_size);
                    linspace(p1_new(2), p2_new(2), discrete_size);
                    linspace(p1_new(3), p2_new(3), discrete_size)];
            end
        end
    end

    % MC update for spatial hashing
    R_hash = numel(P);
    Disc = size(P{1}, 2);

    % Convert 3xD matrices into vectors of length 3D, then stack them.
    pos = reshape( cell2mat(P(:)'), 3*R_hash*Disc, 1 );
    % Build flat pos vector once (3 * total_points, column vector)
    R = num_segments;
    D = discrete_size;
    total_points = R * D;
    %pos = reshape( cell2mat(P(:)') , 3*R*D, 1 );  % 3 x (R*D) -> flatten columnwise

    % Create spatial hash for this time step
    cell_size = 2*proximity_threshold;
    %H = Hash(cell_size, total_points);
    V = VoxelMap(cell_size);
    V = V.create(pos);
    
    %%% ===========================
    %%% MAIN SEGMENT PAIR DETECTION
    %%% ===========================

    %for every segment
    for i = 1:num_segments
        %candidates = [];    % [pi  pj  dist  j_segment]
        candidates = zeros(0, 4);
            
            % === Loop over points in segment i ===
            base_i = (i-1)*D;
            
            for pI = 1:D
                
                global_i = base_i + pI;
                
                % --- Query voxel hash ---
                % --- these are the id's of individual POINTS that are from
                % the discretized segments, however it's an ID for a flat
                % array of all of them
                ids = V.query(pos, global_i, 2*proximity_threshold);
                if isempty(ids), continue; end
                
                % --- Remove self ---
                ids(ids == global_i) = [];
                if isempty(ids), continue; end
                
                % --- Convert global index → (segment jcand, point pj)
                jcand = ceil(double(ids) / double(D));
                mask_j = (jcand > i);   % only consider segments j > i
                
                ids     = ids(mask_j);
                jcand   = jcand(mask_j);
                
                if isempty(ids), continue; end
                
                % --- Skip segments already connected to i ---
                mask_not_connected = ~connected(i, jcand);
                ids     = ids(mask_not_connected);
                jcand   = jcand(mask_not_connected);
                
                if isempty(ids), continue; end
                
                % --- Compute actual distances (vectorized) ---
                pi_vec = pos(3*(global_i-1)+1 : 3*global_i);
                pj_vecs = posreshape(pos, ids);  % helper below
                diffs = pj_vecs - pi_vec;
                dists = sqrt(sum(diffs.^2, 1));
                
                % Filter by threshold
                valid = (dists(:) < proximity_threshold);
                
                %jcand is the fiber # of the j fiber
                if any(valid)

                    pj = double(mod(ids(valid)-1, D) + 1);
                    %takes the point # on segment i for all entries, all j
                    %points that are based on the valid ids, the distances between them,
                    %and segment j it belongs to
                    candidates = [candidates; pI*ones(nnz(valid),1), pj', dists(valid)', jcand(valid)' ];
                end
                
            end
            
            % NOTE ABOUT CANDIDATES
            % The it is semi rare for more than two to ever be in the same
            % vicinity based on the cell size.
            % may add a small if check to see how often that actually
            % happens.

          

            % === Process candidates if any ===
            if ~isempty(candidates)
                % take the closest connecting pair
                % min returns the minimum VALUE and INDEX.
                % we only want index
                %[~, k] = min(candidates(:,3));
                %best = candidates(k,:);
                %  ^^^ ignore ^^^
                
                % Calculate weights (inverse distance, so closer pairs are more likely)
                weights = 1 ./ (candidates(:,3) + 1e-10); % Add small epsilon to avoid division by zero
                weights = weights / sum(weights); % Normalize to probabilities

                % Select a candidate based on weighted probability
                cumulative_weights = cumsum(weights);
                rand_val = rand();
                selected_idx = find(cumulative_weights >= rand_val, 1);

                closest_point_i = candidates(selected_idx, 1);
                closest_point_j = candidates(selected_idx, 2);
                min_dist = candidates(selected_idx, 3);
                segment_j = candidates(selected_idx, 4);
                if isempty(segment_j) || segment_j < 1 || segment_j > num_segments
                    warning('Invalid segment_j detected: %s', mat2str(segment_j));
                    continue;
                end
                
                % Calculate midpoint between the selected points
                % point at segment i, subpoint closet_point_i
                point_i = P{i}(:, closest_point_i);
                point_j = P{segment_j}(:, closest_point_j);
                midpoint = (point_i + point_j) / 2;

                % Move both segments so their connection points meet at the midpoint
                displacement_i = midpoint - point_i;
                displacement_j = midpoint - point_j;

                P{i} = P{i} + displacement_i;
                P{segment_j} = P{segment_j} + displacement_j;

                % record connection
                connected(i, segment_j) = true;
                connected(segment_j, i) = true;
                connections = [connections; i, segment_j, closest_point_i, closest_point_j];
            
         
            
                %Record formation time and discretized points
                junction_formation_time = [junction_formation_time, step * dt];  % Convert step to time
                junction_points = [junction_points; closest_point_i, closest_point_j];


                % Determine connection type for reporting
                %
                %
                is_endpoint_i = (closest_point_i == 1) || (closest_point_i == discrete_size);
                is_endpoint_j = (closest_point_j == 1) || (closest_point_j == discrete_size);

                connection_type = '';
                if is_endpoint_i && is_endpoint_j
                    connection_type = ' (endpoint-endpoint)';
                    new_junc = 1; %1 means end-to-end, or L-junction
                    count_L = count_L + 1;
                    Lform(count_L, 1:2) = 0;
                    % which components came together?
                    if seg_in_comp(i) == 1
                        Lform(count_L, 1) = 1;
                    end
                    if seg_in_comp(segment_j) == 1
                        Lform(count_L, 2) = 1;
                    end
                        
                    %newpart
                    % Calculate interior angle for L-junction
                    if closest_point_i == 1
                        v_i = P{i}(:,2) - P{i}(:,1);  % Point away from endpoint 1
                    else  % closest_point_i == discrete_size
                        v_i = P{i}(:,discrete_size) - P{i}(:,discrete_size-1);  % Point away from last endpoint
                    end

                    if closest_point_j == 1
                        v_j = P{segment_j}(:,2) - P{segment_j}(:,1);
                    else  % closest_point_j == discrete_size
                        v_j = P{segment_j}(:,discrete_size) - P{segment_j}(:,discrete_size-1);
                    end

                    new_angle = acos(dot(v_i, v_j) / (norm(v_i) * norm(v_j)));

                elseif is_endpoint_i || is_endpoint_j
                    connection_type = ' (endpoint-interior)';
                    new_junc = 2; %2 means end-to-middle, or T-junction
                    % For T-junctions, use overall segment directions
                    new_angle = acos(dot((P{i}(:,discrete_size)-P{i}(:,1)),(P{segment_j}(:,discrete_size)-P{segment_j}(:,1)))/(norm(P{i}(:,discrete_size)-P{i}(:,1))*norm(P{segment_j}(:,discrete_size)-P{segment_j}(:,1))));

                else
                    connection_type = ' (interior-interior)';
                    new_junc = 3; %3 means middle-to-middle, or X-junction
                    % For X-junctions, use overall segment directions
                    new_angle = acos(dot((P{i}(:,discrete_size)-P{i}(:,1)),(P{segment_j}(:,discrete_size)-P{segment_j}(:,1)))/(norm(P{i}(:,discrete_size)-P{i}(:,1))*norm(P{segment_j}(:,discrete_size)-P{segment_j}(:,1))));

                end

                theta = [theta new_angle]; %measured in radians

                % New variable: Segment numbers involved in junctions
                seg_in_comp(i) = 1;
                seg_in_comp(segment_j) = 1;

                fprintf('Step %d: Segments %d and %d connected at points %d and %d%s (distance: %.3f -> 0.000) (angle: %.3f)\n', ...
                    step, i, segment_j, closest_point_i, closest_point_j, connection_type, min_dist, new_angle);
                junction_type=[junction_type new_junc];
            end
    end



    

    % Update plots
    if mod(step,plot_interval) == 0 %only plot every "plot_interval" time steps
        step
        countstep=countstep+1;
        for i = 1:num_segments
            set(segment_plots(i), 'XData', P{i}(1,:), ...
                'YData', P{i}(2,:), ...
                'ZData', P{i}(3,:));
        end

        % Visualize connections - draw both connection points
        for conn_idx = 1:size(connections, 1)
            seg1 = connections(conn_idx, 1);
            seg2 = connections(conn_idx, 2);
            point1 = connections(conn_idx, 3);
            point2 = connections(conn_idx, 4);

            % Verify that connection points are still coincident
            pos1 = P{seg1}(:, point1);
            pos2 = P{seg2}(:, point2);
            connection_error = norm(pos1 - pos2);

            if connection_error > 1e-3 % Small tolerance for numerical errors
                fprintf('Warning: Connection %d has error %.6f\n', conn_idx, connection_error);
            end

            % Draw connection points
            plot3(pos1(1), pos1(2), pos1(3), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'red');
            plot3(pos2(1), pos2(2), pos2(3), 'go', 'MarkerSize', 4, 'MarkerFaceColor', 'green');
        end

        drawnow;

        % % Capture frame for storing (optional, if you want to keep this)
        % if countstep == 1
        %     frames = getframe(gcf);
        % else
        %     frames(countstep) = getframe(gcf);
        % end

        % Capture and write frame to video
        frame = getframe(gcf);
        writeVideo(v, frame);

        % Clear connection markers for next frame
        children = get(gca, 'Children');
        markers_to_delete = [];
        for c = 1:length(children)
            if strcmp(get(children(c), 'Marker'), 'o')
                markers_to_delete = [markers_to_delete, c];
            end
        end
        delete(children(markers_to_delete));
    end

end %end of main simulation loop

toc

% Display final statistics
total_connections = size(connections, 1);

fprintf('\nSimulation completed!\n');
fprintf('Total connections formed: %d\n', total_connections);
%fprintf('Frames captured: %d\n', length(frames));

fprintf('number of X-junctions: %d\n',length(find(junction_type==3)));
fprintf('number of T-junctions: %d\n',length(find(junction_type==2)));
fprintf('number of L-junctions: %d\n',length(find(junction_type==1)));

%convert angle from radians to degrees
thetadeg=theta*180/pi;

% %only use the acute angles
% thetadegacute=min(thetadeg, 180 - thetadeg);

%save the data from the run
save -ascii run1_junction_type.dat junction_type
save -ascii run1_junction_anglefull.dat thetadeg
%save -ascii run1_junction_angleacute.dat thetadegacute
save run1_components.mat components
save('run1_parameters.mat', ...
    'L', 'segment_length','num_segments','Tend','dt','proximity_threshold',...
    'discrete_size','D_parallel','D_perp','D_rot');
save -ascii run1_Ljunc.dat Lform
save -ascii run1_nearness.dat nearness
save -ascii run1_junction_time.dat junction_formation_time
save -ascii run1_junction_points.dat junction_points

figure(2)
histogram(junction_type)
title('Histogram of junction types')
savefig('junction_type_run1.fig')


close(v);
disp(['Video saved to: ' videoFileName]);


%visualize histograms of angles
% figure(3)
% histogram(thetadegacute,'BinWidth',5)
% title('Histogram of angles')
% savefig('angles_run1.fig')


% MC update for spatial hashing
function P = posreshape(pos, ids)
    % pos is (3*N)x1
    num = numel(ids);
    P = zeros(3, num);
    for k = 1:num
        idx = ids(k);
        P(:,k) = pos(3*(idx-1)+1 : 3*idx);
    end
end





