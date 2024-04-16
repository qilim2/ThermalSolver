function [full_out,nodes,epoch_keep,Qdot_tot_keep,MESHGRIDS_11,elem_nodes_map] = solver_integrated_R08(bulk_data_file,RK,t_data,plot_switch,rad_switches,env_rad_type,env_data_paths,env_rad_data,s2s_rad_data,workspace_name)
% Other saved values: Qdot_out_keep, Qdot_solar_nodes_keep, Qdot_solar_elems_keep, Qdot_body_nodes_keep, Qdot_body_elems_keep, Qdot_e2e_keep, vS_keep, vBody_keep

% Solver of the transient thermal analysis problem. Takes a bulk data
% file, spits out output temperatures which are both saved and plotted.
% This version sockets with FF to output the temperatures at
% each time step. 

% Version 8.0 completed 4/5/2024

fprintf('---SOLVER RUNNING---\n')

% Importing java libraries to create sockets
import java.net.*
import java.io.*

%% Solver Initialization

% bulk_data_file = path to the bulk data file (string)
% RK = choose RK1 (Forward Euler) or RK4 for stepping
% t_data = holds time inputs for time step size and final time if needed: [dt, tf]
% plot_switch = turn on plot (0 or 1)
% rad_switches = radiation out, environmental radiation, internal radiation vector of 0s (off) and 1s (on)
% env_rad_type = define if the propagator information for the environmental radiation is from a saved file (0) or from a concurrently run FF file (1)
% env_data_paths = path to the saved file or a string vector with: [FreeFlyerPath, PathToMissionPlanFolder, MissionPlanName]
% env_rad_data = input data for environmental radiation (input ~ if not used): [R_body, T_body, albedo_factor]
% s2s_rad_data = input data for internal radiation (input ~ if not used): [reflection_cutoff, ray_count]

if rad_switches(1) == 0
    fprintf('Radiation rejection is turned OFF\n')
elseif rad_switches(1) == 1
    fprintf('Radiation rejection is turned ON\n')
else
    fprintf('Error with rad_switches element 1 input.\n')
    return
end
if rad_switches(2) == 0
    fprintf('Environmental radiation is turned OFF\n');
    t0 = 0; % Initial time
    dt = t_data(1);
    tf = t_data(2);
    t = t0:dt:tf;
    nt = length(t);
    fprintf('Simulation start time: t0 = %f\nSimulation end time: tf = %f\nInitial step size: dt = %f\nNumber of time steps: nt = %i\n',t0,tf,dt,nt);
    vS = 0;
    vBody = 0;
    
elseif rad_switches(2) == 1
    fprintf('Environmental radiation is turned ON\n')
    R_body = env_rad_data(1);
    T_body = env_rad_data(2);
    albedo_factor = env_rad_data(3);

    AF = albedo_factor; % Albedo factor (recommended to use 0.36 per Thornton book)

    if env_rad_type == 0
        input_file = env_data_paths;
        % Read input file
        [epoch_vec,vS_mat,vBody_mat,eclipse_vec] = input_file_parser_R01(input_file);
        t = epoch_vec - epoch_vec(1);
        t = t(2:end); % Vector of times
        dt_vec = epoch_vec(2:end) - epoch_vec(1:end-1); % Vector of time steps
        nt = length(dt_vec);
        fprintf('Set number of time steps: %i\n',nt)
        fprintf('Data successfully read in from text file.\n')
    elseif env_rad_type == 1
        FreeFlyerPath = env_data_paths{1};
        PathToMissionPlanFolder = env_data_paths{2};
        MissionPlanName = env_data_paths{3};
        t = 0;
    else
        fprintf('Error with rad_env_type in solver_integrated input.\n');
        return
    end
else
    fprintf('Error with rad_switches element 2 input.\n');
    return
end
if rad_switches(3) == 0
    fprintf('Surface-to-surface heat exchange is turned OFF\n');
elseif rad_switches(3) == 1
    fprintf('Surface-to-surface heat exchange is turned ON\n');
    reflection_cutoff = s2s_rad_data(1);
    ray_count = s2s_rad_data(2);

    cutoff = reflection_cutoff; % Cutoff value for reflections
else
    fprintf('Error with rad_switches element 3 input.\n');
end


SBC = 5.670374419*10^(-8); % Stefan-Boltzmann Constant [W/m2/K4]



%% Bulk Data Processing
% Read bulk data file
[overview, MESHGRIDS_11, SURFACES_11, ~, ~, THERMOOPTICAL_11, ELEMENTS_11, CONDUCTANCES_11, HEATLOADS_11] = bulk_data_parser_R02(bulk_data_file);

% Overview
node_count = overview(1);
elem_count = overview(2);

% Prepare mass and temp
nodes = MESHGRIDS_11(:,1);
Cth = MESHGRIDS_11(:,2);
T = MESHGRIDS_11(:,3);

% Rename for clarity or brevity
G = CONDUCTANCES_11;
Qdot_HL = HEATLOADS_11;

% Store temperatures
% T_keep1 = Inf*ones(length(nodes),nt);
T_keep1(:,1) = T;
MESHGRIDS_12 = MESHGRIDS_11; % Hold temperatures and nodal information for the current step

%%% Prep surface values
surf_gtop = SURFACES_11(:,10);
surf_gtop_sn = SURFACES_11(:,12);
surf_gbot = SURFACES_11(:,13);
surf_gbot_sn = SURFACES_11(:,15);

%%% Prep element info
elem_ID = (1:elem_count)';
elem_nodes = ELEMENTS_11(:,2:5);
[~,elem_nodes_map] = ismember(elem_nodes,nodes); % Map nodes within elements
elem_areas = ELEMENTS_11(:,6);
elem_normals = ELEMENTS_11(:,7:9);
elem_alpha = zeros(height(ELEMENTS_11),2);
elem_epsilon = zeros(height(ELEMENTS_11),2);
elem_gtop = zeros(height(ELEMENTS_11),1); elem_gtop_role = elem_gtop; elem_gtop_sn = elem_gtop;
elem_gbot = zeros(height(ELEMENTS_11),1); elem_gbot_role = elem_gbot; elem_gbot_sn = elem_gbot;
for ii = 1:elem_count
    surf_index = find(ELEMENTS_11(ii,1)==SURFACES_11(:,1));
    elem_alpha(ii,1) = THERMOOPTICAL_11(SURFACES_11(surf_index,4),2); % Topside
    elem_alpha(ii,2) = THERMOOPTICAL_11(SURFACES_11(surf_index,5),2); % Bottomside
    elem_epsilon(ii,1) = THERMOOPTICAL_11(SURFACES_11(surf_index,4),3); % Topside
    elem_epsilon(ii,2) = THERMOOPTICAL_11(SURFACES_11(surf_index,5),3); % Bottomside
    elem_gtop(ii,1) = SURFACES_11(surf_index,10); % Group ID of topside for this element
    elem_gtop_role(ii,1) = SURFACES_11(surf_index,11); % Role (0=ignored, 1=active+blocker, 2=blocker) of topside for this element
    elem_gtop_sn(ii,1) = SURFACES_11(surf_index,12); % Boolean indicator of spacenode visibility for the topside for this element
    elem_gbot(ii,1) = SURFACES_11(surf_index,13); % Group ID of bottomside for this element
    elem_gbot_role(ii,1) = SURFACES_11(surf_index,14); % Role (0=ignored, 1=active+blocker, 2=blocker) of bottomside for this element
    elem_gbot_sn(ii,1) = SURFACES_11(surf_index,15); % Boolean indicator of spacenode visibility for the bottomside for this element
end
elem_gtop_sn = logical(elem_gtop_sn); % Convert to boolean/logical
elem_gbot_sn = logical(elem_gbot_sn); % Convert to boolean/logical

%% Radiation Preparation
if sum(rad_switches) > 0
    %%% Element center coordinates
    elem_center = zeros(elem_count,3);
    for ii = 1:elem_count
        elem_center(ii,:) = rect_elem_center_R01(ii,ELEMENTS_11,MESHGRIDS_11);
    end
    
    %%% Prep for rad_out
    % Isolate and map elements to go into rad_out function
    elem_top_rad_out = NaN; elem_bot_rad_out = NaN; % Initialize values
    pre_rad_out_top = NaN; pre_rad_out_bot = NaN; % Initialize values
    if sum(elem_gtop_role==1) ~= 0
        elem_top_rad_out = elem_ID(elem_gtop_role==1); % Get IDs of only active elements. Needed for later remapping.
        % pre_rad_out_top = SBC*elem_epsilon(elem_gtop_role==1,1).*elem_areas(elem_gtop_role==1); % Prepare topside coefficient.
        % pre_rad_out_top = repelem(pre_rad_out_top/4,4); % pre_rad_out_coefficients is by elements. The values in this vector need to be repeated four times each (and divided by 4) to match the number of nodes.
        % elem_top_rad_out_mapped_nodes = elem_nodes_map(elem_gtop_role==1,:); % Get relevant nodes with new indices
        % elem_top_rad_out = elem_gtop_role==1;
        pre_rad_out_top = SBC*elem_epsilon(:,1).*elem_areas; % Elemental
        pre_rad_out_top = pre_rad_out_top.*(elem_gtop_role==1); % Set values for elems that don't radiate to zero
        pre_rad_out_top = repelem(pre_rad_out_top/4,4); % Convert to nodal

    end
    if sum(elem_gbot_role==1) ~= 0
        elem_bot_rad_out = elem_ID(elem_gbot_role==1); % Get IDs of only active elements. Needed for later remapping.
        % pre_rad_out_bot = SBC*elem_epsilon(elem_gbot_role==1,2).*elem_areas(elem_gbot_role==1); % Prepare bottomside_coefficient
        % pre_rad_out_bot = repelem(pre_rad_out_bot/4,4); % pre_rad_out_coefficients is by elements. The values in this vector need to be repeated four times each (and divided by 4) to match the number of nodes.
        % elem_bot_rad_out_mapped_nodes = elem_nodes_map(elem_gbot_role==1,:); % Get relevant nodes with new indices
        % elem_bot_rad_out = elem_gbot_role==1;
        pre_rad_out_bot = SBC*elem_epsilon(:,2).*elem_areas;
        pre_rad_out_bot = pre_rad_out_bot.*(elem_gbot_role==1);
        pre_rad_out_bot = repelem(pre_rad_out_bot/4,4);
    end
    
    
    %%% Prep for rad_environmental
    % Environmental radiation should implement radiation groups as there should be
    % shadowing of surfaces. Shadowing should be calculated on an element
    % basis if computation speed allows.
    % Environmental radiation is assumed ON if the space node boolean is activated
    % for the surface. This means that the surface, even if it is a blocker or
    % ignored within a radiation group can experience incoming external
    % radiation. The group role only determines the object's role within the
    % radiation group. Therefore, all elements part of a group within which at
    % least one elem/surf experiences external radiation need to be passed into
    % the external radiation function.
    if rad_switches(2) == 1 || rad_switches(3) == 1
        % Find groups with at least one surface that sees sn
        radgroups_rad_env = unique([surf_gtop(surf_gtop_sn == 1);surf_gbot(surf_gbot_sn == 1)]); % Get list of rad group IDs
        g = cell(1,length(radgroups_rad_env));
        for ii = 1:length(radgroups_rad_env) % Iterate through each group
            % Get current group ID
            gID = radgroups_rad_env(ii);
        
            % Get list of elements
            f1 = find(elem_gtop==gID);
            f2 = find(elem_gbot==gID);
            elems = unique([f1;f2]);
        
            % Store info
            h.ID = gID;
            h.elems = elems;
            h.tsn = boolean(elem_gtop_sn(elems));
            h.bsn = boolean(elem_gbot_sn(elems));
            h.trole = elem_gtop_role(elems);
            h.brole = elem_gbot_role(elems);
            h.centers = elem_center(elems,:);
            
            % Store in group cell array
            g{ii} = h;
        end
        
        % Prepare rotation matrices for each element to transform the plane and any
        % points on the plane to a 2D projection (goes from s/c body frame to
        % element frame. To go from element frame to s/c body frame left multiply by 
        % the transpose of this matrix.)
        elem_rotms = cell(1,elem_count);
        for ii = 1:elem_count
            unit_norm = [0,0,1]; % Create a vector to align with
            if sum(cross(elem_normals(ii,:),unit_norm),"all")==0
                if dot(unit_norm,elem_normals(ii,:)) >= 0
                    R = eye(3);
                else
                    R = -eye(3);
                end
            else
                R = create_rotmat_R01(elem_normals(ii,:),unit_norm);
            end
            elem_rotms{ii} = R;
        end
        
        % Get coordinates of each point in the element
        ecoords = cell(1, elem_count); % Preallocate cell array
        for ii = 1:elem_count
            current_nodes = elem_nodes_map(ii,:);
            ecoords{ii} = [MESHGRIDS_11(current_nodes(1),4:6)',MESHGRIDS_11(current_nodes(2),4:6)',...
                MESHGRIDS_11(current_nodes(3),4:6)',MESHGRIDS_11(current_nodes(4),4:6)'];
        end
    end
    
    %%% Prep for radiation_internal
    if rad_switches(3) == 1
        % Radiation proportion calculations
        fprintf('Calculating internal view factors and RPM.\n');
        RPM = cell(1,length(g));
        for ii = 1:length(g)
            RPM{ii} = radiation_internal_prep_R01(g{ii},ecoords,ray_count,elem_rotms,elem_epsilon,cutoff);
        end
        assignin(workspace_name,'RPMs',RPM)
    end
end
fprintf('Solver initialization complete.\n');

% check = sum(RPM{1},2)

%% FF connection
if rad_switches(2) == 1 && env_rad_type == 1
    
    showGUI = 0; % Open FF GUI if desired
    terminationCode = '10191980'; % Tells Matlab when to end Java read
    
    %%% Create 2 sockets per FF instance
    % Try to open socket 1 using the OS-determined port
    try
        socketServer(1) = ServerSocket(0); % JAVA command
    catch err
        error(strcat('Error: Unable to open socket. ', getReport(err)))
    end
    % Get the port number used for the above socket
    portNum(1) = socketServer(1).getLocalPort(); % JAVA command
    
    % Try to open socket 2 using the OS-determined port
    try
        socketServer(2) = ServerSocket(0); % JAVA command
    catch err
        error(strcat('Error: Unable to open socket. ', getReport(err)))
    end
    % Get the port number used for the above socket
    portNum(2) = socketServer(2).getLocalPort(); % JAVA command 
    
    %%% Launch FF Instance    
    % Create command-line string to be executed
    if boolean(showGUI)
       commandString = strcat('"', FreeFlyerPath, 'FreeFlyer.exe"', ...
           ' -r -mp "', PathToMissionPlanFolder, MissionPlanName,'"', ...
           sprintf(' -ui %d -ui %d -ui %s &', portNum(1), portNum(2), terminationCode));
    else
       commandString = strcat('"', FreeFlyerPath, 'FF.exe"', ...
           ' -mp "', PathToMissionPlanFolder, MissionPlanName,'"', ...
           sprintf(' -ui %d -ui %d -ui %s  &', portNum(1), portNum(2), terminationCode));
    end
     
    % Execute command at command-line.  Other execution commands are dos() and
    % !.  However, both have subtle differences from system()
    disp(commandString)
    system( commandString );
    
    % Wait for the FF client instance to connect.
    socketServer(1).setSoTimeout(int16(10000)); % JAVA command
    
    % Open first socket.
    try
        socketClient(1) = socketServer(1).accept(); % JAVA command
    catch err
        error('Error: Unable to accept MissionPlan as client')
    end
    
    socketServer(2).setSoTimeout(int16(10000)); % JAVA command
    
    %open second socket
    try
        socketClient(2) = socketServer(2).accept(); % JAVA command
    catch err
        error('Error: Unable to accept MissionPlan as client')
    end
    
    %%% Create write/read buffer
    outputStream = DataOutputStream(socketClient(1).getOutputStream()); % JAVA command
    inputStream  = DataInputStream(socketClient(2).getInputStream()); % JAVA command
    inputReader = BufferedReader(InputStreamReader(inputStream)); % JAVA command

    %%% Get initial set of data
    % Get information from orbital solver
    % Request data from FF client
    orbit_data_request = 1;
    outputStream.writeDouble(orbit_data_request); % JAVA command
    
    %%%  Wait for FF Instances to finish running     
    %{
    This very simply listens to the data being input through the stream and
    only increments when non-null data is read in.  This reader will only break
    when it receives the end-loop command of '10191980'.
    %}
    index = 1;
    while 1  
        JAVAdata = inputReader.readLine(); % JAVA command
        rawData{ index } = char(JAVAdata);    
        
        if ( strcmp(terminationCode,rawData{ index }) == 1 )
            break
        end
        
        index = index + 1;
        
    end
    
    % Convert data
    received_data = str2num(rawData{1});
    epoch_keep = received_data(1);
end



%% Sim
fprintf('Beginning simulation...\n');
% fprintf('Total Time [s] = %2.1f\nCurrent time [s] = ',t(end)) % Set up timer output
t0 = 0;
time_out = fprintf('%2.1f',t0); % Display initial time
loop_number = 2;
end_loop = 0;
while end_loop ~= 1

    Tn1 = T_keep1(:,loop_number-1); % Current temperature of all nodes
    
    %%% Get information from orbital solver
    if rad_switches(2) == 1 && env_rad_type == 0
        vS = vS_mat(loop_number-1,:);
        vBody = vBody_mat(loop_number-1,:);
        eclipse = eclipse_vec(loop_number-1,:);
        dt = dt_vec(loop_number-1);
    elseif rad_switches(2) == 1 && env_rad_type == 1
        % Request data from FF client
        orbit_data_request = 1;
        outputStream.writeDouble(orbit_data_request); % JAVA command

        %%%  Wait for FF Instances to finish running     
        %{
        This very simply listens to the data being input through the stream and
        only increments when non-null data is read in.  This reader will only break
        when it receives the end-loop command of '10191980'.
        %}
        index = 1;
        while 1  
            JAVAdata = inputReader.readLine(); % JAVA command
            rawData{ index } = char(JAVAdata);    
            
            if ( strcmp(terminationCode,rawData{ index }) == 1 )
                break
            end
            
            index = index + 1;
            
        end
        
        % Convert data
        received_data = str2num(rawData{1});
        epoch_keep(loop_number) = received_data(1); % Keep epoch
        t(loop_number) = received_data(1) - epoch_keep(1); % Get sim time
        dt = epoch_keep(loop_number) - epoch_keep(loop_number-1);
        vS = received_data(2:4)*1000;
        vBody = received_data(5:7)*1000;
        eclipse = received_data(8);
        end_loop = received_data(9);
        if end_loop == 1
            break
        end
    end
    vS_keep(loop_number-1,:) = vS;
    vBody_keep(loop_number-1,:) = vBody;


    %%% Radiation out
    % Calculate Qdot_out for each node given the immediate area around it.
    % Use pre_rad_out_top and bot and divide to each node in element.
    % This requires finding each node per element and doing an assignment.
    % Check nodal masses function for calculating the area per node.
    % Sum total radiation out and store for later checks.
    if rad_switches(1) == 1 || rad_switches(3) == 1
        % Input active elements into rad_out function
        Qdot_out_nodes_top = zeros(node_count,1); Qdot_out_nodes_bot = zeros(node_count,1); % Initialize both in case one is not solved.
        Qdot_out_elems_top = zeros(elem_count,1); Qdot_out_elems_bot = zeros(elem_count,1); 
        if ~isnan(elem_top_rad_out(1)) % Check if NaN by checking first value only.
            % Work with element topsides
            % Qdot_out_top = radiation_out_R01(Tn1,pre_rad_out_top,elem_top_rad_out_mapped_nodes,node_count); % Note: this should match new indices of nodes. Output is per node.
            % [Qdot_out_nodes_top,Qdot_out_elems_top] = radiation_out_R02(Tn1,pre_rad_out_top,elem_top_rad_out_mapped_nodes,node_count,elem_top_rad_out,elem_count);
            % [Qdot_out_nodes_top,Qdot_out_elems_top] = radiation_out_R02(Tn1,pre_rad_out_top,elem_nodes_map,node_count,elem_top_rad_out,elem_count);
            [Qdot_out_nodes_top,Qdot_out_elems_top] = radiation_out_R03(Tn1,pre_rad_out_top,elem_nodes_map);
        end
        if ~isnan(elem_bot_rad_out(1)) % Check if NaN by checking first value only.
            % Work with element bottomsides
            % Qdot_out_bot = radiation_out_R01(Tn1,pre_rad_out_bot,elem_bot_rad_out_mapped_nodes,node_count); % Note: this should match new indices of nodes. Output is per node.
            % [Qdot_out_nodes_bot,Qdot_out_elems_bot] = radiation_out_R02(Tn1,pre_rad_out_bot,elem_bot_rad_out_mapped_nodes,node_count,elem_bot_rad_out,elem_count);
            [Qdot_out_nodes_bot,Qdot_out_elems_bot] = radiation_out_R03(Tn1,pre_rad_out_bot,elem_nodes_map);
        end
        Qdot_out_nodes = Qdot_out_nodes_top + Qdot_out_nodes_bot; % Sum
        % Qdot_out_elems = Qdot_out_elems_top + Qdot_out_elems_bot;
    else
        Qdot_out_nodes = zeros(node_count,1);
    end
    
    Qdot_out_keep(:,loop_number-1) = Qdot_out_nodes;

    %%% Environmental radiation
    % Calculate Qdot_in for each element given its view of the space node.
    % Input vectors to sun and earth. Check against top and bottom sides of
    % elements that are part of active surfaces.
    if rad_switches(2) == 1
        % Vector to sun, vector to orbital body
        vSrep = repmat(vS,elem_count,1); % Prepare by repeating the vS vector for each element to make elem_count x 3 matrix
        vbodyrep = repmat(vBody,elem_count,1); % Prepare by repeating the vbody vector for each element to make elem_count x 3 matrix
        
        % Solar incidence angles
        solar_incidence_angles = vector_angles_R03(elem_normals,vSrep); % Get angles of incidence for each element
        solar_angle_bool = (solar_incidence_angles < pi/2); % Get boolean indicating if seeing top or bottomside
        elem_solar_proj_area = abs(elem_areas.*cos(solar_incidence_angles)); % Get the projected area relative to the angle
        
        % Orbital body incidence angles and view factors (topside)
        theta_t = vector_angles_R03(elem_normals,vbodyrep); % Get angles of incidence for each element
        sintheta_t = sin(theta_t);
        costheta_t = cos(theta_t);
        cottheta_t = cot(theta_t);
    
        % Orbital body incidence angles and view factors (bottomside)
        theta_b = vector_angles_R03(-elem_normals,vbodyrep);
        sintheta_b = sin(theta_b);
        costheta_b = cos(theta_b);
        cottheta_b = cot(theta_b);
    
        eta = asin(R_body/norm(vBody)); % Maximum angle subtended by Earth (angle between nadir vector of s/c and line drawn from satellite to horizon of Earth)
        H = norm(vBody)/R_body;
        % elem_body_proj_area_t = abs(elem_areas.*costheta_t); % Get the projected area relative to the topside angle
        % elem_body_proj_area_b = abs(elem_areas.*costheta_b); % Get the projected area relative to the bottomside angle
        
        % Create boolean vectors indicating which sides are visible
        illum_topside_only = theta_t < (pi/2 - eta);
        illum_bottomside_only = theta_t > (pi/2 + eta);
        illum_both_sides = theta_t <= pi/2 + eta & theta_t > pi/2 - eta;
        illum_t = illum_topside_only | illum_both_sides; % Boolean true if topside is illuminated either way
        illum_b = illum_bottomside_only | illum_both_sides;
    
        % Compose view factor vectors
        elem_body_t_VF = zeros(size(illum_topside_only));
        simple_t_VF = costheta_t/(H^2);
        complex_t_VF = (2/pi)*(pi/4 - 0.5*asin(sqrt(H^2-1)./(H*sintheta_t))+(1/(2*H^2)).*(costheta_t.*acos(-sqrt(H^2-1).*cottheta_t)-sqrt(H^2-1).*sqrt(1-H^2*(costheta_t.^2))));
        elem_body_t_VF(illum_topside_only) = simple_t_VF(illum_topside_only);
        elem_body_t_VF(illum_both_sides) = complex_t_VF(illum_both_sides);
    
        elem_body_b_VF = zeros(size(illum_bottomside_only));
        simple_b_VF = costheta_b/(H^2);
        complex_b_VF = (2/pi)*(pi/4 - 0.5*asin(sqrt(H^2-1)./(H*sintheta_b))+(1/(2*H^2)).*(costheta_b.*acos(-sqrt(H^2-1).*cottheta_b)-sqrt(H^2-1).*sqrt(1-H^2*(costheta_b).^2)));
        elem_body_b_VF(illum_bottomside_only) = simple_b_VF(illum_bottomside_only);
        elem_body_b_VF(illum_both_sides) = complex_b_VF(illum_both_sides);
    
        Qdot_solar_elems = zeros(elem_count,length(radgroups_rad_env));
        Qdot_albedo_elems = zeros(elem_count,length(radgroups_rad_env));
        Qdot_IR_elems = zeros(elem_count,length(radgroups_rad_env));
        parfor ii = 1:length(radgroups_rad_env)
            [Qdot_solar_elems(:,ii),Qdot_albedo_elems(:,ii),Qdot_IR_elems(:,ii),S_to_keep(ii)] = radiation_environmental_R03(T_body,AF,vS,vBody,g{ii},elem_normals,ecoords,elem_rotms,elem_alpha,elem_epsilon,solar_angle_bool,elem_solar_proj_area,illum_t,illum_b,elem_body_t_VF,elem_body_b_VF,eclipse,elem_areas,R_body);
        end
        Qdot_solar_elems = sum(Qdot_solar_elems,2);
        Qdot_albedo_elems = sum(Qdot_albedo_elems,2);
        Qdot_IR_elems = sum(Qdot_IR_elems,2);
        S_keep(loop_number-1,1) = S_to_keep(1);
    else
        Qdot_solar_elems = zeros(elem_count,1);
        Qdot_albedo_elems = zeros(elem_count,1);
        Qdot_IR_elems = zeros(elem_count,1);
    end
    Qdot_solar_nodes = elem2nodes_R01(Qdot_solar_elems,elem_nodes_map,node_count);
    Qdot_solar_nodes_keep(:,loop_number-1) = Qdot_solar_nodes;
    Qdot_solar_elems_keep(:,loop_number-1) = Qdot_solar_elems;
    Qdot_albedo_nodes = elem2nodes_R01(Qdot_albedo_elems,elem_nodes_map,node_count);
    Qdot_albedo_nodes_keep(:,loop_number-1) = Qdot_albedo_nodes;
    Qdot_albedo_elems_keep(:,loop_number-1) = Qdot_albedo_elems;
    Qdot_IR_nodes = elem2nodes_R01(Qdot_IR_elems,elem_nodes_map,node_count);
    Qdot_IR_nodes_keep(:,loop_number-1) = Qdot_IR_nodes;
    Qdot_IR_elems_keep(:,loop_number-1) = Qdot_IR_elems;


    %%% Internal radiation
    if rad_switches(3) == 1
        % Note: Should be entirely balanced within the satellite. Include a check. May need to check against Qdot_out.
        Qdot_e2e = zeros(elem_count,length(radgroups_rad_env));
        parfor ii = 1:length(radgroups_rad_env)
            Qdot_e2e(:,ii) = radiation_internal_R01(g{ii},RPM{ii},elem_count,Qdot_out_elems_top,Qdot_out_elems_bot);
        end
        Qdot_e2e = sum(Qdot_e2e,2);
    else
        Qdot_e2e = zeros(elem_count,1);
    end
    Qdot_e2e_nodes = elem2nodes_R01(Qdot_e2e,elem_nodes_map,node_count);
    Qdot_e2e_keep(:,loop_number-1) = Qdot_e2e_nodes;

    %%% Sum heats on each element (include env rad results) and separate out to nodes
    % Note: already separated out to nodes for radiation out.
    % Qdot_elems = Qdot_env + Qdot_e2e;
    % Qdot_in = elem2nodes_R01(Qdot_elems,elem_nodes_map,node_count);

    %%% Sum nodal heating rates
    if rad_switches(1) == 0
        Qdot_out_nodes = 0;
    end
    if rad_switches(2) == 0
        Qdot_solar_nodes = 0;
        Qdot_body_nodes = 0;        
    end
    if rad_switches(3) == 0
        Qdot_e2e_nodes = 0;
    end
    Qdot_tot = Qdot_solar_nodes + Qdot_albedo_nodes + Qdot_IR_nodes + Qdot_e2e_nodes - Qdot_out_nodes;
    Qdot_tot_keep(:,loop_number-1) = Qdot_tot;

    
    %%% Solver
    if RK == 1
        % Heat transfer with RK1
        deltaT = Tn1.'-Tn1;
        k1 = (Qdot_HL + Qdot_tot + sum(G.*deltaT,2))./Cth;
        T_change = k1*dt;
        T = Tn1 + T_change;
        MESHGRIDS_12(:,3) = T;
        T_keep1(:,loop_number) = T;
    elseif RK == 4
        % Heat transfer with RK4
        deltaT = Tn1.'-Tn1;
        k1 = (Qdot_HL + Qdot_tot + sum(G.*deltaT,2))./Cth;
        deltaT = (Tn1+k1/2*dt).'-(Tn1+k1/2*dt);
        k2 = (Qdot_HL + Qdot_tot + sum(G.*deltaT,2))./Cth;
        deltaT = (Tn1+k2/2*dt).'-(Tn1+k2/2*dt);
        k3 = (Qdot_HL + Qdot_tot + sum(G.*deltaT,2))./Cth;
        deltaT = (Tn1+k3*dt).'-(Tn1+k3*dt);
        k4 = (Qdot_HL + Qdot_tot + sum(G.*deltaT,2))./Cth;
        T_change = (k1+2*k2+2*k3+k4)*dt/6;
        T = Tn1 + T_change;
        MESHGRIDS_12(:,3) = T;
        T_keep1(:,loop_number) = T;
    elseif RK == 15
        fun = @(T,Tn) (Qdot_HL + Qdot_tot + sum(G.*(Tn.'-Tn),2))./Cth;
        tspan = [t(loop_number-1),t(loop_number)];
        [~,T] = ode15s(fun,tspan,Tn1);
        MESHGRIDS_12(:,3) = T(end,:)';
        T_keep1(:,loop_number) = T(end,:)';
    elseif RK == 45
        fun = @(T,Tn) (Qdot_HL + Qdot_tot + sum(G.*(Tn.'-Tn),2))./Cth;
        tspan = [t(loop_number-1),t(loop_number)];
        [~,T] = ode45(fun,tspan,Tn1);
        MESHGRIDS_12(:,3) = T(end,:)';
        T_keep1(:,loop_number) = T(end,:)';
    elseif RK == 89
        fun = @(T,Tn) (Qdot_HL + Qdot_tot + sum(G.*(Tn.'-Tn),2))./Cth;
        tspan = [t(loop_number-1),t(loop_number)];
        [~,T] = ode89(fun,tspan,Tn1);
        MESHGRIDS_12(:,3) = T(end,:)';
        T_keep1(:,loop_number) = T(end,:)';
    else
        fprintf('ERROR (solver_integrated): Invalid input for RK. Use either RK=1, RK=4, or RK=45.\n')
        return
    end

    %%% Plot
    if plot_switch == 1
        figure(3)
        surface_post_visualization_R02(MESHGRIDS_12,ELEMENTS_11) % Plot
    end
    
    %%% Gif
    %exportgraphics(gca,"checkout_surface_01_01.gif","Append",true)
    
    %%% Timer
    if rad_switches(2) == 1 && env_rad_type == 1
        fprintf(repmat('\b',1,time_out))
        time_out = fprintf('%2.1f',epoch_keep(loop_number)-epoch_keep(1));
    elseif rad_switches(2) == 1 && env_rad_type == 0
        fprintf(repmat('\b',1,time_out))
        time_out = fprintf('%2.1f',t(loop_number));
        epoch_keep = epoch_vec;
    else
        fprintf(repmat('\b',1,time_out))
        time_out = fprintf('%2.1f',t(loop_number));
        epoch_keep = t;
    end

    %%% Update loop number and check for end
    % Note that end_loop will handled by FF if using env_rad_type == 1
    if rad_switches(2) == 0
        % fprintf('Test 1. nt = %i, loop number = %i\n',nt,loop_number);
        if loop_number >= nt
            end_loop = 1;
        end
    elseif rad_switches(2) == 1
        if env_rad_type == 0
           if loop_number >= nt
               end_loop = 1;
           end

           % The ELSE is handled by FF passing end_loop=1 to matlab
        end
    end
    
    loop_number = loop_number + 1;
    % disp(end_loop)
end

fprintf('\n')

%% Output
full_out = T_keep1;
if rad_switches(1) == 1
    assignin(workspace_name,'Qdot_out_keep',Qdot_out_keep);
end
if rad_switches(2) == 1 
    assignin(workspace_name,'Qdot_solar_nodes_keep',Qdot_solar_nodes_keep);
    assignin(workspace_name,'Qdot_solar_elems_keep',Qdot_solar_elems_keep);
    assignin(workspace_name,'Qdot_albedo_nodes_keep',Qdot_albedo_nodes_keep);
    assignin(workspace_name,'Qdot_albedo_elems_keep',Qdot_albedo_elems_keep);
    assignin(workspace_name,'Qdot_IR_nodes_keep',Qdot_IR_nodes_keep);
    assignin(workspace_name,'Qdot_IR_elems_keep',Qdot_IR_elems_keep);
    if env_rad_type == 1
        epoch_keep(end) = [];
        assignin(workspace_name,'S_keep',S_keep);
        assignin(workspace_name,'vS_keep',vS_keep);
        assignin(workspace_name,'vBody_keep',vBody_keep);
    elseif env_rad_type == 0
        epoch_keep(end) = [];
        assignin(workspace_name,'S_keep',S_keep);
        assignin(workspace_name,'vS_keep',vS_keep);
        assignin(workspace_name,'vBody_keep',vBody_keep);
    end
end
if rad_switches(3) == 1
    assignin(workspace_name,'Qdot_e2e_keep',Qdot_e2e_keep);
end
end