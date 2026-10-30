function kinetic_net_non(example, epsilon, tau, t_end, nx_inner, nv)
%KINETIC_NET_NON  Reference solver for the kinetic BGK model on a network.
%
%   KINETIC_NET_NON(example, epsilon, tau, t_end, nx_inner, nv, save_data)
%   solves the scaled kinetic BGK equation for the perturbation f of the
%   distribution F = M + epsilon*f on every edge of a network,
%
%       df/dt + v/epsilon * df/dx = ( (M[F] - M)/epsilon - f ) / (epsilon^2 * tau),
%
%   with kinetic coupling conditions at the nodes. It serves as a reference
%   for the diffusion-limit solvers (global_netnondc_hf, global_netnondc_m).
%
%   NOTE: This solver is implemented in a deliberately naive way (explicit
%   upwind transport with dt ~ epsilon*dx, full velocity grid) to obtain
%   accurate reference values. It is therefore slow for small epsilon.
%   A considerably faster approximation will be added to this repository
%   (work in progress).
%
%   Inputs:
%     example   - ID of the network topology (see setup_network_topology)
%     epsilon   - scaling parameter (Knudsen number)
%     tau       - relaxation parameter
%     t_end     - final time
%     nx_inner  - number of inner grid cells per edge (edges are [0,1])
%     nv        - number of velocity cells on [-8, 8] (default: 500)
%
%   Output:
%     The moments rhoK, uK, pK, TK are saved to
%     Data/non_kinex<...>.mat. During the run, a checkpoint is written to
%     the same file every 10 minutes (wall-clock time).
%
%   Dependencies: setup_network_topology.m, couplinghf.m

arguments
    example (1,1) double
    epsilon (1, 1) double
    tau (1, 1) double
    t_end (1, 1) double
    nx_inner (1, 1) double
    nv = 500
end

%% Initialization & Setup
[nE, nN, N, pB, uB, rhoB] = setup_network_topology(example);
% N{iN}.iE   : edges attached to node iN
% N{iN}.sign : -1 if the edge ends at the node (x = 1), +1 if it starts there (x = 0)
% N{iN}.ad   : 0 for exterior (boundary) edges, > 0 for inner edges

% Kinetic transmission matrix at the nodes (outgoing = C3 * incoming):
% example 0: two edges, full transmission; otherwise: three edges, equal split
if example == 0
    C3 = ones(2) - eye(2);
else
    C3 = (ones(3) - eye(3)) / 2;
end

% Spatial discretization (one ghost cell at each edge end)
x_L = 0;
x_R = 1;
nbc = 1;
nx = nx_inner + 2*nbc;
dx = (x_R - x_L) / nx_inner;
inner = nbc+1 : nx_inner+nbc;

% Velocity grid (symmetric around 0) and Maxwellian
v_min = -8;
v_max = 8;
dv = (v_max - v_min) / nv;
v = v_min + dv/2 : dv : v_max - dv/2;
M = (1/sqrt(2*pi)) * exp(-v.^2 / 2);

vm = v < 0;
vp = v > 0;
v_vp = v(vp);
v_vm = v(vm);
num_vm = sum(vm);

%% Initial and boundary states (one value per edge)
p0 = ones(1, nE);
u0 = zeros(1, nE);
rho0 = ones(1, nE);

rho_bc_L_mac = rho0;
u_bc_L_mac   = u0;
p_bc_L_mac   = p0;

rho_bc_R_mac = rhoB;
u_bc_R_mac   = uB;
p_bc_R_mac   = pB;

% Assign macroscopic boundary values depending on edge orientation
for iN = 1:nN
    for ie = 1:length(N{iN}.iE)
        IE = N{iN}.iE(ie);
        aD = N{iN}.ad(ie);

        if N{iN}.sign(ie) < 0
            rho_bc_R_mac(IE) = rho0(IE);
            u_bc_R_mac(IE)   = u0(IE);
            p_bc_R_mac(IE)   = p0(IE);

            if aD <= 0
                rho_bc_L_mac(IE) = rhoB(IE);
                u_bc_L_mac(IE)   = uB(IE);
                p_bc_L_mac(IE)   = pB(IE);
            end
        end
    end
end

% Linearized kinetic boundary distributions (rows: edges, columns: velocities)
BC_L = (rho_bc_L_mac' + u_bc_L_mac'.*v + 0.5*(v.^2 - 1).*(p_bc_L_mac' - rho_bc_L_mac')) .* M;
BC_R = (rho_bc_R_mac' + u_bc_R_mac'.*v + 0.5*(v.^2 - 1).*(p_bc_R_mac' - rho_bc_R_mac')) .* M;

% Half-range moments (0th, 1st, 2nd) of the distribution entering each node
rho1 = zeros(1, nE);
u1   = zeros(1, nE);
p1   = zeros(1, nE);

for iN = 1:nN
    for ie = 1:length(N{iN}.iE)
        IE = N{iN}.iE(ie);
        if N{iN}.sign(ie) < 0
            rho1(IE) = sum(BC_L(IE, vp)) * dv;
            u1(IE)   = sum(BC_L(IE, vp) .* v_vp) * dv;
            p1(IE)   = sum(BC_L(IE, vp) .* (v_vp.^2)) * dv;
        else
            rho1(IE) = sum(BC_R(IE, vm)) * dv;
            u1(IE)   = sum(BC_R(IE, vm) .* v_vm) * dv;
            p1(IE)   = sum(BC_R(IE, vm) .* (v_vm.^2)) * dv;
        end
    end
end

% Coupled macroscopic state, used for a well-prepared initial condition
g = [rho1; u1; p1];
[p_coup, u_coup] = couplinghf(example, g, nE, nN, N);

%% Well-prepared initial distribution f(x, v, edge)
f = zeros(nx, nv, nE);
macro_rho_init = (sum(BC_L, 2) + sum(BC_R, 2)) * dv / 2;

for iE = 1:nE
    well_prep_init = (macro_rho_init(iE) + v.*u_coup(iE) + 0.5*(v.^2 - 1).*(p_coup(iE) - macro_rho_init(iE))) .* M;

    f(:, :, iE) = ones(nx, 1) * well_prep_init;
    f(1, :, iE) = BC_L(iE, :);
    f(end, :, iE) = BC_R(iE, :);
end

%% Kinetic solver loop
time_kin = 0;

dataname = "Data/non_kinex" + example + "t" + t_end + "eps" + epsilon + "tau" + tau + "nx" + nx_inner + "nv" + nv + ".mat";
temp_dataname = "Data/tempnon_kinex" + example + "t" + t_end + "eps" + epsilon + "tau" + tau + "nx" + nx_inner + "nv" + nv + ".mat";

save_interval_seconds = 600;   % checkpoint interval (wall-clock seconds)
last_save_time = tic;

rhoK = zeros(nx, nE);
uK   = zeros(nx, nE);
TK   = zeros(nx, nE);
pK   = zeros(nx, nE);

while (time_kin < t_end)
    % CFL condition of the explicit transport step
    dt_kin = epsilon * dx / max(abs(v_min), abs(v_max));
    dt_kin = min(dt_kin, t_end - time_kin);

    dt_eps_dx = dt_kin / (dx * epsilon);
    relax_rate1 = dt_kin / (epsilon^3 * tau);
    relax_rate2 = 1 + dt_kin / (epsilon^2 * tau);

    % 1. Boundary conditions: incoming velocities from the BCs, outgoing extrapolated
    for iE = 1:nE
        f(1:nbc, :, iE)  = f(nbc+1, :, iE);
        f(1:nbc, vp, iE) = BC_L(iE, vp);
        f(end, :, iE)    = f(end-nbc, :, iE);
        f(end, vm, iE)   = BC_R(iE, vm);
    end

    % 2. Kinetic node coupling: distribute the distributions leaving each
    %    edge into the node onto the other edges (fliplr maps v -> -v
    %    on the symmetric velocity grid)
    for iN = 1:nN
        num_edges = length(N{iN}.iE);
        fin = zeros(num_edges, num_vm);

        for iE_idx = 1:num_edges
            IE = N{iN}.iE(iE_idx);
            if N{iN}.sign(iE_idx) > 0
                fin(iE_idx, :) = f(nbc+1, vm, IE);
            else
                fin(iE_idx, :) = fliplr(f(end-1, vp, IE));
            end
        end

        fout = C3(1:num_edges, 1:num_edges) * fin;

        for iE_idx = 1:num_edges
            IE = N{iN}.iE(iE_idx);
            if N{iN}.sign(iE_idx) > 0
                f(1, vp, IE) = fliplr(fout(iE_idx, :));
                f(1, vm, IE) = f(2, vm, IE);
            else
                f(end, vm, IE) = fout(iE_idx, :);
                f(end, vp, IE) = f(end-1, vp, IE);
            end
        end
    end

    % 3. Transport (first-order upwind) and BGK relaxation
    for iE = 1:nE
        f(inner, vp, iE) = f(inner, vp, iE) - dt_eps_dx * (f(inner, vp, iE) - f(inner-1, vp, iE)) .* v(vp);
        f(inner, vm, iE) = f(inner, vm, iE) - dt_eps_dx * (f(inner+1, vm, iE) - f(inner, vm, iE)) .* v(vm);

        % Moments of the perturbation f
        f_in = f(inner, :, iE);
        rhoK(inner,iE) = sum(f_in, 2) * dv;
        uK(inner, iE)   = sum(f_in .* v, 2) * dv;
        TK(inner, iE)   = sum(f_in .* (v.^2 - 1), 2) * dv;

        % Density, velocity and temperature of the full distribution F = M + epsilon*f
        rhoF = 1 + epsilon * rhoK(inner,iE);
        uF   = (epsilon * uK(inner,iE)) ./ rhoF;
        TF   = 1 + epsilon * TK(inner,iE) ./ (1 + epsilon * rhoK(inner,iE)) - epsilon^2 * (uK(inner,iE).^2) ./ (1 + epsilon * rhoK(inner,iE)).^2;

        % Local Maxwellian of F
        M_F = (rhoF ./ sqrt(2 * pi * TF)) .* exp(-0.5 * (v - uF).^2 ./ TF);

        % Relaxation step, implicit in f
        f(inner, :, iE) = (f_in + relax_rate1 * (M_F - M)) / relax_rate2;
    end

    % Periodic checkpoint (written to a temporary file first, then moved)
    if toc(last_save_time) > save_interval_seconds

        fprintf('Saving checkpoint at simulation time: %f...\n', time_kin);

        for iE = 1:nE
            rhoK(:, iE) = sum(f(:,:,iE), 2) * dv;
            uK(:, iE)   = sum(f(:,:,iE) .* v, 2) * dv;
            TK(:, iE)   = sum(f(:,:,iE) .* (v.^2 - 1), 2) * dv;
            pK(:, iE)   = sum(f(:,:,iE) .* v.^2, 2) * dv;
        end

        save(temp_dataname, "example", "epsilon", "tau", "time_kin", "t_end", "nx_inner", "nv", "rhoK", "uK", "pK", "TK");
        movefile(temp_dataname, dataname, 'f');

        last_save_time = tic;
    end

    time_kin = time_kin + dt_kin;
    %fprintf('Time: %f\n', time_kin);
end


%% Final moments and data export
rhoK = zeros(nx, nE);
uK   = zeros(nx, nE);
TK   = zeros(nx, nE);
pK   = zeros(nx, nE);

for iE = 1:nE
    rhoK(:, iE) = sum(f(:,:,iE), 2) * dv;
    uK(:, iE)   = sum(f(:,:,iE) .* v, 2) * dv;
    TK(:, iE)   = sum(f(:,:,iE) .* (v.^2 - 1), 2) * dv;
    pK(:, iE)   = sum(f(:,:,iE) .* v.^2, 2) * dv;
end

save(dataname, "example", "epsilon", "tau", "t_end", "nx_inner", "nv", "rhoK", "uK", "pK", "TK");

end