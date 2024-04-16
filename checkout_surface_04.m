close all; clear; clc;
format short
tic

%% Pre-processing plot switch
pre_plotting = 1;

%% Initialization
MESHGRIDS_1 = []; % Global holding matrix of all nodes in mesh
surface_count_current = 1; % Current index of surfaces (counts up per number of surfaces, used only as indexing)
SURFACES_1 = {}; % Global holding cell array of surface properties
SURFACES_2 = {}; % Global holding cell array of nodes owned by surfaces and surface nodal definitions
THERMOPHYSICAL_1 = ["thermophysical name","density (rho)","conductivity (k)","specific heat (cp)"]; % Global holding matrix of thermophysical properties
THERMOOPTICAL_1 = ["thermo-optical name","absorptivity (alpha)","emissivity (epsilon)"]; % Global holding matrix of thermooptical properties
CONDUCTANCES_1 = []; % Global holding matrix of conductance values.
CONDUCTANCES_2 = {}; % Cell array that holds contactor information.

%% Surface mesh (create)
% Surface mesh inputs
surface_number = 1;
surf_method = 1; % 0 = global, 1 = relative to origin
origin = [-0.15;-0.05;0.10];
vec_x11 = [0.30;0;0];
vec_x12 = [0;0.10;0];
n_x1 = 7;
n_x2 = 3;
start_node = 1;
[MESHGRIDS_1,SURFACES_2] = surface_mesh_create_R01(surface_number,origin,surf_method,vec_x11,vec_x12,n_x1,n_x2,start_node,MESHGRIDS_1,SURFACES_2);

surface_number = 2;
surf_method = 1;
origin = [-0.15;-0.05;-0.10];
vec_x21 = [0.30;0;0];
vec_x22 = [0;0;0.20];
n_x1 = 7;
n_x2 = 5;
start_node = 22;
[MESHGRIDS_1,SURFACES_2] = surface_mesh_create_R01(surface_number,origin,surf_method,vec_x21,vec_x22,n_x1,n_x2,start_node,MESHGRIDS_1,SURFACES_2);

surface_number = 3;
surf_method = 1;
origin = [0.15;-0.05;-0.10];
vec_x31 = [0;0.10;0];
vec_x32 = [0;0;0.20];
n_x1 = 3;
n_x2 = 5;
start_node = 57;
[MESHGRIDS_1,SURFACES_2] = surface_mesh_create_R01(surface_number,origin,surf_method,vec_x31,vec_x32,n_x1,n_x2,start_node,MESHGRIDS_1,SURFACES_2);

surface_number = 4;
surf_method = 1;
origin = [0.15;0.05;-0.10];
vec_x41 = [-0.30;0;0];
vec_x42 = [0;0;0.20];
n_x1 = 7;
n_x2 = 5;
start_node = 72;
[MESHGRIDS_1,SURFACES_2] = surface_mesh_create_R01(surface_number,origin,surf_method,vec_x41,vec_x42,n_x1,n_x2,start_node,MESHGRIDS_1,SURFACES_2);

surface_number = 5;
surf_method = 1;
origin = [-0.15;0.05;-0.10];
vec_x51 = [0;-0.10;0];
vec_x52 = [0;0;0.20];
n_x1 = 3;
n_x2 = 5;
start_node = 107;
[MESHGRIDS_1,SURFACES_2] = surface_mesh_create_R01(surface_number,origin,surf_method,vec_x51,vec_x52,n_x1,n_x2,start_node,MESHGRIDS_1,SURFACES_2);

surface_number = 6;
surf_method = 1;
origin = [-0.15;0.05;-0.10];
vec_x61 = [0.30;0;0];
vec_x62 = [0;-0.10;0];
n_x1 = 7;
n_x2 = 3;
start_node = 122;
[MESHGRIDS_1,SURFACES_2] = surface_mesh_create_R01(surface_number,origin,surf_method,vec_x61,vec_x62,n_x1,n_x2,start_node,MESHGRIDS_1,SURFACES_2);

fprintf('Surface mesh created.\n')

%% Thermophysical Properties
thermophysical_name = "aluminum";
rho = 2711; % Density [kg/m3]
k = 179.690; % Thermal conductivity [W/m/K]
% k = 237;
cp = 896; % Specific heat [J/kg/K]

THERMOPHYSICAL_1 = thermophysical_create_edit_R01(thermophysical_name,rho,k,cp,THERMOPHYSICAL_1);
fprintf('Thermophsyical properties created.\n')

%% Thermo-optical Properties
thermooptical_name_01 = "Thermo_Optical_01";
alpha = 0.15; % Solar absorptivity
epsilon = 0.8; % Infrared emissivity
THERMOOPTICAL_1 = thermooptical_create_edit_R01(thermooptical_name_01,alpha,epsilon,THERMOOPTICAL_1);

thermooptical_name_02 = "Thermo_Optical_02";
alpha = 0.3; % Solar absorptivity
epsilon = 0.5; % Infrared emissivity
THERMOOPTICAL_1 = thermooptical_create_edit_R01(thermooptical_name_02,alpha,epsilon,THERMOOPTICAL_1);

fprintf('Thermo-optical properties created.\n')

%% Surface Properties
T_init = 273.15; % Initial temperature [K]
thickness = 0.003; % Surface thickness [m]

nodes2edit = -1;
surfaces2edit = [1,2,3,4,5,6];

% Assign surface properties to nodes
[MESHGRIDS_1,SURFACES_1] = surface_properties_assign_R02(1,"aluminum","Thermo_Optical_01","Thermo_Optical_01",vec_x11,vec_x12,thickness,MESHGRIDS_1,SURFACES_1,SURFACES_2); 
[MESHGRIDS_1,SURFACES_1] = surface_properties_assign_R02(2,"aluminum","Thermo_Optical_01","Thermo_Optical_01",vec_x21,vec_x22,thickness,MESHGRIDS_1,SURFACES_1,SURFACES_2); 
[MESHGRIDS_1,SURFACES_1] = surface_properties_assign_R02(3,"aluminum","Thermo_Optical_01","Thermo_Optical_01",vec_x31,vec_x32,thickness,MESHGRIDS_1,SURFACES_1,SURFACES_2); 
[MESHGRIDS_1,SURFACES_1] = surface_properties_assign_R02(4,"aluminum","Thermo_Optical_01","Thermo_Optical_01",vec_x41,vec_x42,thickness,MESHGRIDS_1,SURFACES_1,SURFACES_2); 
[MESHGRIDS_1,SURFACES_1] = surface_properties_assign_R02(5,"aluminum","Thermo_Optical_01","Thermo_Optical_01",vec_x51,vec_x52,thickness,MESHGRIDS_1,SURFACES_1,SURFACES_2); 
[MESHGRIDS_1,SURFACES_1] = surface_properties_assign_R02(6,"aluminum","Thermo_Optical_01","Thermo_Optical_01",vec_x61,vec_x62,thickness,MESHGRIDS_1,SURFACES_1,SURFACES_2); 
fprintf('Surface properties assigned to nodes.\n')

% Assign thermal masses to nodes
MESHGRIDS_1 = set_nodal_thermal_mass_R01(surfaces2edit,thermophysical_name,THERMOPHYSICAL_1,MESHGRIDS_1,SURFACES_1,SURFACES_2);
fprintf('Thermal masses assigned to nodes.\n')

% Assign initial temperature to nodes
MESHGRIDS_1 = set_nodal_temperature_R02(T_init,nodes2edit,surfaces2edit,MESHGRIDS_1,SURFACES_2);
fprintf('Thermal masses assigned to nodes.\n')

%% Create HL
% Apply sample heat load along edge of surface 1
HL_nodes = 1:21;
HL_01 = zeros(height(MESHGRIDS_1),1);
HL = 0;
HL_node_count = length(HL_nodes);
HL_01(HL_nodes) = HL/HL_node_count;

fprintf('Heat load(s) created.\n')

%% Define Elements
ELEMENTS_1 = elements_define_R01(MESHGRIDS_1,SURFACES_2);
fprintf('Elements defined.\n')

%% Nodal Lengths
SURFACELENGTHS_1 = Inf*ones(height(MESHGRIDS_1)); % Initialize holding matrix after MESHGRIDS_1 is created
SURFACELENGTHS_1 = nodal_lengths_R02(surfaces2edit,MESHGRIDS_1,ELEMENTS_1,SURFACELENGTHS_1);
fprintf('Surface Lengths calculated.\n')

%% Surface Conductances
CONDUCTANCES_1 = surface_conductance_R02(1,CONDUCTANCES_1,SURFACES_1,SURFACES_2,THERMOPHYSICAL_1,SURFACELENGTHS_1);
CONDUCTANCES_1 = surface_conductance_R02(2,CONDUCTANCES_1,SURFACES_1,SURFACES_2,THERMOPHYSICAL_1,SURFACELENGTHS_1);
CONDUCTANCES_1 = surface_conductance_R02(3,CONDUCTANCES_1,SURFACES_1,SURFACES_2,THERMOPHYSICAL_1,SURFACELENGTHS_1);
CONDUCTANCES_1 = surface_conductance_R02(4,CONDUCTANCES_1,SURFACES_1,SURFACES_2,THERMOPHYSICAL_1,SURFACELENGTHS_1);
CONDUCTANCES_1 = surface_conductance_R02(5,CONDUCTANCES_1,SURFACES_1,SURFACES_2,THERMOPHYSICAL_1,SURFACELENGTHS_1);
CONDUCTANCES_1 = surface_conductance_R02(6,CONDUCTANCES_1,SURFACES_1,SURFACES_2,THERMOPHYSICAL_1,SURFACELENGTHS_1);
fprintf('Conductance matrix created.\n')

%% Merge nodes or contact
merge = 1;
if merge == 1 % Merge nodes at an edge
    % Check mass beforehand
    mass_check = sum(MESHGRIDS_1(:,2));
    fprintf('Thermal mass before merge: %f\n',mass_check)

    % Do merge
    range_tol = 1e-4;
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(1,31,2,33,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s1 South to s2 North
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(1,32,3,33,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s1 East to s3 North
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(1,33,4,33,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s1 North to s4 North
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(1,34,5,33,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s1 West to s5 North
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(2,32,3,34,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s2 East to s3 West
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(3,32,4,34,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s3 East to s4 West
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(4,32,5,34,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s4 East to s5 West
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(5,32,2,34,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s5 East to s2 West
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(2,31,6,33,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s2 South to s6 North
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(3,31,6,32,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s3 South to s6 East
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(4,31,6,31,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s4 South to s6 South
    [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(5,31,6,34,range_tol,0,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2); % s5 South to s6 West

    % Check mass again
    mass_check = sum(MESHGRIDS_1(:,2));
    fprintf('Thermal mass after merge: %f\n',mass_check)

elseif merge == 2 % Create contactor at the edge instead
    priority_switch = 1;
    CONDUCTANCES_2 = contactor_create_R01(1,1,30,2,30,100,3,0,0.5,priority_switch,MESHGRIDS_1,SURFACES_1,SURFACES_2,CONDUCTANCES_2);
end

%% Rad Groups
% Choose which surfaces and what sides experience external radiation (into satellite)
% (ID,topbot,group_ID,role,sn,SURFACES_1)
% ID - Surface ID
% topbot - Select both (1), topside (2), bottomside (3)
% group_ID - Group number
% role - Role of the surface (0=ignored, 1=active+blocker, 2=blocker)
% sn - Boolean for turning on/off spacenode vision (0=off, 1=on)
SURFACES_1 = rad_groups_assign_R01(1,2,1,1,1,SURFACES_1); % Assign topside of surface 1 to external group. Can see sn.
SURFACES_1 = rad_groups_assign_R01(1,3,2,1,0,SURFACES_1); % Assign bottomside of surface 1 to internal group. Cannot see sn.
SURFACES_1 = rad_groups_assign_R01(2,2,1,1,1,SURFACES_1); % Assign topside of surface 2 to external group. Can see sn.
SURFACES_1 = rad_groups_assign_R01(2,3,2,1,0,SURFACES_1); % Assign bottomside of surface 2 to internal group Cannot see sn.
SURFACES_1 = rad_groups_assign_R01(3,2,1,1,1,SURFACES_1); % Assign topside of surface 1 to external group. Can see sn.
SURFACES_1 = rad_groups_assign_R01(3,3,2,1,0,SURFACES_1); % Assign bottomside of surface 1 to internal group. Cannot see sn.
SURFACES_1 = rad_groups_assign_R01(4,2,1,1,1,SURFACES_1); % Assign topside of surface 2 to external group. Can see sn.
SURFACES_1 = rad_groups_assign_R01(4,3,2,1,0,SURFACES_1); % Assign bottomside of surface 2 to internal group Cannot see sn.
SURFACES_1 = rad_groups_assign_R01(5,2,1,1,1,SURFACES_1); % Assign topside of surface 1 to external group. Can see sn.
SURFACES_1 = rad_groups_assign_R01(5,3,2,1,0,SURFACES_1); % Assign bottomside of surface 1 to internal group. Cannot see sn.
SURFACES_1 = rad_groups_assign_R01(6,2,1,1,1,SURFACES_1); % Assign topside of surface 2 to external group. Can see sn.
SURFACES_1 = rad_groups_assign_R01(6,3,2,1,0,SURFACES_1); % Assign bottomside of surface 2 to internal group Cannot see sn.

fprintf('Radiation groups defined.\n');

%% Plotting
% Plot nodes
marker_area = 500;
back_switch = 1;
normal_switch = 1;
surfaces2plot = [1,2,3,4,5,6];

if pre_plotting == 1
    figure(1)
    surface_pre_visualization_R01(surfaces2plot,back_switch,marker_area,MESHGRIDS_1,SURFACES_2)
    xlabel('X'); ylabel('Y');zlabel('Z');
    
    figure(2)
    mesh_visualization_R01(back_switch,normal_switch,MESHGRIDS_1,ELEMENTS_1);
    xlabel('X'); ylabel('Y');zlabel('Z');
    
    fprintf('Pre-visualization complete.\n')
    
end

%% Reset timer
toc
tic

%% Save
bulk_data='checkout_04.blk';
HEATLOADS_1 = HL_01;

bulk_data_process_R02(bulk_data,MESHGRIDS_1,SURFACES_1,SURFACES_2,THERMOPHYSICAL_1,THERMOOPTICAL_1,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2,HEATLOADS_1)

%% Solve
% clear

bulk_data_file='checkout_04.blk'; % Name/path of the bulk data file
RK = 89; % 1-Forward Euler, 4-RK4, 15-ode15s, 45-ode45, 89-ode89
dt = 30; % Time step size [s]
tf = 32400; % Final time step
t_data = [dt,tf]; % Inputs for time step and final time if not defined by the input file
plot_switch = 0; % Turn off/on plotting during the simulation
rad_switches = [1,1,1]; % Turn on/off radiation out, environmental radiation, surface-to-surface radiation
env_rad_type = 1; % Switch between saved data (0) or concurrent FF file (1)
input_file='C:\Users\qilim2\Box\Thesis\FreeFlyer Missions\output_test_06.txt'; % Name/path of saved data
FreeFlyerPath = 'C:\Program Files\a.i. solutions, Inc\FreeFlyer 7.8.0.56841 (64-Bit)\'; % Path to the FF executable
PathToMissionPlanFolder = 'C:\Users\qilim2\Box\Thesis\FreeFlyer Missions\'; % Path to the folder where the mission plan is saved
MissionPlanName = 'checkout_interface_06.MissionPlan'; % Name of the mission plan file
if env_rad_type == 0
    env_data_paths = input_file;
elseif env_rad_type == 1
    env_data_paths = {FreeFlyerPath,PathToMissionPlanFolder,MissionPlanName};
end
R_body = 6378e3; % Radius of body [m]
T_body = 255; % Temperature of body [k]
AF = 0.36; % Albedo factor. Use 0.36 per Thornton book
env_rad_data = [R_body,T_body,AF];
cutoff = 1e-3; % Cutoff value for reflections
ray_count = 2e3; % Temporary value for ray count
s2s_rad_data = [cutoff,ray_count];
workspace_name = 'base';

[full_out,node_IDs,epoch_keep,Qdot_tot_keep,MESHGRIDS_11,elem_nodes_map] = solver_integrated_R08(bulk_data_file,RK,t_data,plot_switch,rad_switches,env_rad_type,env_data_paths,env_rad_data,s2s_rad_data,workspace_name);
fprintf('Total number of timesteps: %i\n',length(epoch_keep));
t_keep = epoch_keep-epoch_keep(1);

%% End of main script
toc

%% Plot
figure(4)
title('Nodal temperature over time')
plot(t_keep,full_out(11,:));
grid on; hold on;
plot(t_keep,full_out(39,:));
plot(t_keep,full_out(88,:));
% plot(epoch_keep-epoch_keep(1),full_out(1,:));
% plot(epoch_keep-epoch_keep(1),full_out(28,:));
legend('11 (+z)','39 (-y)','132 (*88) (-z)');
xlim([t_keep(1),t_keep(end)]);
xlabel('Time (s)'); ylabel('Temperature (K)');

figure(5)
title('Solar flux over time (element on -z face)')
plot(t_keep(1:188),Qdot_solar_elems_keep(85,1:188)/0.0025);
grid on;
xlabel('Time (s)'); ylabel('Flux (W/m2)');
xlim([t_keep(1),t_keep(188)]);

figure(6)
title('Albedo flux over time (element on +z face)')
plot(t_keep(1:188),Qdot_albedo_elems_keep(9,1:188)/0.0025);
grid on; 
xlabel('Time (s)'); ylabel('Flux (W/m2)');
xlim([t_keep(1),t_keep(188)]);

figure(7)
title('IR flux over time (element on +z face)')
plot(t_keep(1:188),Qdot_IR_elems_keep(9,1:188)/0.0025);
xlabel('Time (s)'); ylabel('Flux (W/m2)');
xlim([t_keep(1),t_keep(188)]);

figure(8)
timestep_number = 1079;
surface_post_visualization_R03(timestep_number,MESHGRIDS_11,full_out,elem_nodes_map)


%% Save
save('results_04_03.mat')

