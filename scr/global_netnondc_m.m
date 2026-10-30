function global_netnondc_m(example, tau, t_end, nx_inner, nv)
%GLOBAL_NETNONDC_M  Diffusion limit of a linear kinetic model on a network
%   with Maxwellian boundary conditions at the nodes.
%
%   GLOBAL_NETNONDC_M(example, tau, t_end, nx_inner, nv, save_data) solves
%   the convection-diffusion equation for z = p - a^2*rho on every edge of a
%   network. Edges are coupled at the nodes by optimized Robin transmission
%   conditions and solved with a Schwarz domain-decomposition iteration.
%   Coupling data at the nodes is obtained from the first and second
%   half-range velocity moments (Maxwellian node conditions, see couplingm).
%
%   Inputs:
%     example   - ID of the network topology (see setup_network_topology)
%     tau       - diffusion coefficient
%     t_end     - final time
%     nx_inner  - number of inner grid cells per edge (edges are [0,1])
%     nv        - number of velocity cells on [-8, 8] (default: 500)
%
%   Output:
%     The esults are saved to Data/min_non_mex<...>.mat.
%
%   Dependencies: setup_network_topology.m, couplingm.m

arguments
    example (1,1) double
    tau (1, 1) double
    t_end (1, 1) double
    nx_inner (1, 1) double
    nv (1,1) double = 500
end

%% Initialization & Setup
delta_2 = (3* sqrt(2* pi))/(2* (pi-2));   % constant in the Robin transmission conditions
a = sqrt(3);                               % sound speed
[nE, nN, N, pB, uB, rhoB] = setup_network_topology(example);
% N{iN}.iE   : edges attached to node iN
% N{iN}.sign : -1 if the edge ends at the node (x = 1), +1 if it starts there (x = 0)
% N{iN}.ad   : 0 for exterior (boundary) edges, > 0 for inner edges

% Spatial discretization (one ghost cell at each edge end)
x_L = 0;
x_R = 1;
nbc = 1;
nx = nx_inner + 2*nbc;
dx = (x_R - x_L) / nx_inner;
inner = nbc+1 : nx_inner+nbc;

% Velocity grid and Maxwellian
v_min = -8;
v_max = 8;
dv = (v_max - v_min) / nv;
v = v_min + dv/2 : dv : v_max - dv/2;
M = (1/sqrt(2*pi)) * exp(-v.^2 / 2);

%% Initial and boundary states (one value per edge)
p0 = ones(1, nE);
u0 = zeros(1, nE);
rho0 = ones(1, nE);

rho_bc_L_mac = rho0; u_bc_L_mac = u0; p_bc_L_mac = p0;
rho_bc_R_mac = rhoB; u_bc_R_mac = uB; p_bc_R_mac = pB;

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

% Linearized kinetic distributions (rows: edges, columns: velocities)
BC_R = (rho_bc_R_mac' + v.*u_bc_R_mac' + 0.5*(v.^2 - 1).*(p_bc_R_mac' - rho_bc_R_mac')) .* M;
BC_L = (rho_bc_L_mac' + v.*u_bc_L_mac' + 0.5*(v.^2 - 1).*(p_bc_L_mac' - rho_bc_L_mac')) .* M;
initial = (rho0' + v.*u0' + 0.5*(v.^2 - 1).*(p0' - rho0')) .* M;

vm = v < 0; vp = v > 0;
v_vm = v(vm); v_vp = v(vp);

% Half-range moments (1st, 2nd) of the distribution entering each node
u1 = zeros(1, nE); p1 = zeros(1, nE);

for iN = 1:nN
    for ie = 1:length(N{iN}.iE)
        IE = N{iN}.iE(ie);
        if N{iN}.sign(ie) < 0
            u1(IE)   = sum(BC_L(IE, vp) .* v_vp) * dv;
            p1(IE)   = sum(BC_L(IE, vp) .* (v_vp.^2)) * dv;
        else
            u1(IE)   = sum(BC_R(IE, vm) .* v_vm) * dv;
            p1(IE)   = sum(BC_R(IE, vm) .* (v_vm.^2)) * dv;
        end
    end
end

%% Initial density (inner cells from the initial state, ghost cells from the BCs)
rhoI = zeros(nx, nE);
rhoI_initial = sum(initial, 2)' * dv;
rhoI_BC_L = sum(BC_L, 2)' * dv;
rhoI_BC_R = sum(BC_R, 2)' * dv;

for iE = 1:nE
    rhoI(:, iE)   = rhoI_initial(iE);
    rhoI(1, iE)   = rhoI_BC_L(iE);
    rhoI(end, iE) = rhoI_BC_R(iE);
end

% Maxwellian node coupling -> constant pressure and velocity per edge
g = [u1; p1];
[p_coup, u_coup, rho_coup] = couplingm(example, g, nE, nN, N);

pI = repmat(p_coup, nx, 1);
uI = repmat(u_coup, nx, 1);
z = pI - a^2 * rhoI;

% Dirichlet values for z at the exterior edge ends
for iN = 1:nN
    for ie = 1:length(N{iN}.iE)
        IE = N{iN}.iE(ie);
        aD = N{iN}.ad(ie);
        if aD == 0
            if N{iN}.sign(ie) < 0
                z(1, IE) = pI(1, IE) - a^2 * rho_coup(IE);
            elseif N{iN}.sign(ie) > 0
                z(end, IE) = pI(end, IE) - a^2 * rho_coup(IE);
            end
        end
    end
end

%% Transmission factors of the discrete scheme
time_incomp = 0;
dt_incomp = dx;
root = zeros(2, nE);        % roots of the characteristic equation per edge
kappa = zeros(nE, 2);       % transmission factors (col 1: x = 0, col 2: x = 1)
alpha_opt = zeros(nE, nN);  % Robin parameters (edge, node)

for iN= 1:nN
    for ie = 1:length(N{iN}.iE)
        IE = N{iN}.iE(ie);
        u_val = uI(1, IE);
        term_sqrt = sqrt(dx^2 + 4*tau*dt_incomp + dt_incomp^2 * u_val^2);

        root(1, IE) = ((dx^2/dt_incomp + 2*tau) - (dx/dt_incomp) * term_sqrt) / (2*tau - u_val*dx);
        root(2, IE) = ((dx^2/dt_incomp + 2*tau) + (dx/dt_incomp) * term_sqrt) / (2*tau - u_val*dx);

        % Exterior edges: factor determined by the Dirichlet condition
        if N{iN}.ad(ie) == 0
            factor = (root(2, IE) * root(1, IE)^nx - root(1, IE) * root(2, IE)^nx) / (root(1, IE)^nx - root(2, IE)^nx);
            if N{iN}.sign(ie) < 0
                kappa(IE, 2) = factor;
            else
                kappa(IE, 1) = factor;
            end
        end
    end
end

% Optimize alpha node by node (one parameter per attached edge end) instead
% of all parameters at once. The transmission factor of an inner edge depends
% on the Robin condition at its other node, so kappa is updated after each
% node optimization (Gauss-Seidel sweep over the nodes) until alpha converges.
alpha_opt = ones(nE, nN);              % start value 1 (0 gives a tiny fminsearch simplex)
for iN = 1:nN
    kappa = update_kappa(kappa, alpha_opt, iN, root, uI, tau, dx, nx, N, delta_2);
end
options = optimset('TolX', 1e-10, 'TolFun', 1e-8, 'MaxFunEvals', 5000, 'MaxIter', 5000);
for sweep = 1:20
    alpha_prev = alpha_opt;  rho_nodes = zeros(1, nN);
    for iN = 1:nN
        objective = @(al) node_spectral_radius(al, tau, dx, kappa, uI, iN, N, delta_2);
        [al, rho_nodes(iN)] = fminsearch(objective, alpha_opt(N{iN}.iE, iN), options);
        alpha_opt(N{iN}.iE, iN) = al;
        kappa = update_kappa(kappa, alpha_opt, iN, root, uI, tau, dx, nx, N, delta_2);
    end
    dalpha = max(abs(alpha_opt - alpha_prev), [], 'all');
    fprintf('alpha sweep %d: max. nodal spectral radius %.2e, max. change in alpha %.2e\n', ...
            sweep, max(rho_nodes), dalpha);
    if dalpha < 1e-6, break; end
end


%% System matrices (implicit Euler, central differences), factorized once
zeta = tau / dx^2;
tri_cells = cell(nE, 1);
tri_dec   = cell(nE, 1);

for IE = 1:nE
    u_val = uI(1, IE) / (dx*2);

    main_diag = (1/dt_incomp + 2 * zeta) * ones(nx_inner, 1);
    sub_diag  = (-zeta - u_val) * ones(nx_inner, 1);
    sup_diag  = (-zeta + u_val) * ones(nx_inner, 1);

    T = speye(nx);   % ghost rows default to identity (Dirichlet)
    T(inner, inner) = spdiags([sub_diag, main_diag, sup_diag], -1:1, nx_inner, nx_inner);
    T(1+nbc, 1)     = -zeta - u_val;
    T(end-nbc, end) = -zeta + u_val;
    tri_cells{IE}   = T;
end

% Replace ghost rows at node-attached ends by the Robin conditions
for iN = 1:nN
    for ie = 1:length(N{iN}.iE)
        IE = N{iN}.iE(ie);
        u_val = uI(1, IE) / 2;
        term_alpha = alpha_opt(IE, iN) * delta_2 / 2;
        if N{iN}.sign(ie) < 0
            tri_cells{IE}(end, end)   = term_alpha + u_val - tau/dx;
            tri_cells{IE}(end, end-1) = tau/dx + term_alpha + u_val;
        else
            tri_cells{IE}(1, 1) = term_alpha - u_val - tau/dx;
            tri_cells{IE}(1, 2) = tau/dx + term_alpha - u_val;
        end
    end
end

for IE = 1:nE
    tri_dec{IE} = decomposition(tri_cells{IE});
end

% Constant terms of the nodal coupling conditions
C = zeros(nN, 1);
Cij = zeros(3, 3, nN);
for iN = 1:nN
    sign_uI = N{iN}.sign .* uI(1, N{iN}.iE);
    pI_val  = pI(1, N{iN}.iE);

    C(iN, 1) = 2 * pI_val * sign_uI';
    Cij(:,:,iN) = 3 * (sign_uI' - sign_uI) + delta_2 * (pI_val' - pI_val);
end

%% Time loop
dt_incomp = dx;
for n = 1:round(t_end/dt_incomp)
    i = 0;
    z_old = z;

    % Schwarz iteration: until convergence, at most 100 iterations
    dz = inf;
    while i < 100 && dz > 1e-12 * max(1, max(abs(z(:))))
        z_prev = z;
        rhs = zeros(nx, nE);
        for iN = 1:nN
            for ie = 1:length(N{iN}.iE)
                IE = N{iN}.iE(ie);
                rhs(inner, IE) = z_old(inner, IE) / dt_incomp;
                index = setdiff(1:length(N{iN}.iE), ie);   % neighbouring edges at this node

                % Keep Dirichlet value at the exterior end
                if N{iN}.ad(ie) == 0
                    if N{iN}.sign(ie) < 0
                        rhs(1, IE) = z(1, IE);
                    else
                        rhs(end, IE) = z(end, IE);
                    end
                end

                % Robin data from the neighbouring edges (previous iterate)
                sign_idx = N{iN}.sign(index);
                iE_idx   = N{iN}.iE(index);

                col_idx0 = 1 + (sign_idx < 0) * (nx - 1);   % ghost cell at the node
                col_idx1 = 2 + (sign_idx < 0) * (nx - 3);   % adjacent inner cell

                lin0 = col_idx0 + (iE_idx - 1) * nx;
                lin1 = col_idx1 + (iE_idx - 1) * nx;

                s0 = tau/dx + alpha_opt(IE, iN)*delta_2/4 + sign_idx .* uI(1, iE_idx)/2;
                s1 = -tau/dx + alpha_opt(IE, iN)*delta_2/4 + sign_idx .* uI(1, iE_idx)/2;

                boundary_term = s0 * z(lin0)' + s1 * z(lin1)' + alpha_opt(IE, iN) * sum(Cij(ie, :, iN))/2 + C(iN, 1);

                if N{iN}.sign(ie) < 0
                    rhs(end, IE) = boundary_term;
                else
                    rhs(1, IE) = boundary_term;
                end
            end
        end

        % Independent solves on all edges
        for IE = 1:nE
            z(:, IE) = tri_dec{IE} \ rhs(:, IE);
        end
        dz = max(abs(z(:) - z_prev(:)));
        i = i + 1;
    end
    if dz > 1e-8
        warning('Schwarz iteration did not converge after %d iterations (change %.1e).', i, dz);
    end

    % Recover density from z
    rhoI = (pI - z) / a^2;
    time_incomp = time_incomp + dt_incomp;

    %fprintf('Time: %f\n', time_incomp);
end

% Save results
dataname = "Data/min_non_mex" + example + "t" + t_end +"tau" + tau+ "nx" + nx_inner+ "nv" + nv + ".mat";
save(dataname, "example", "tau", "t_end", "nx_inner", "nv", "rhoI", "uI", "pI", "z", "rhoB", "uB", "pB");
end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function rho = node_spectral_radius(al, tau, dx, kappa, uI, iN, N, delta_2)
%NODE_SPECTRAL_RADIUS  Spectral radius of the Robin iteration at node iN
%   for the parameters al (one value per attached edge end).
n = numel(N{iN}.iE);  Mi = zeros(n);
for ie = 1:n
    IE = N{iN}.iE(ie);  s = N{iN}.sign(ie);
    idx = setdiff(1:n, ie);
    kj = kappa(sub2ind(size(kappa), N{iN}.iE(idx), 1 + (N{iN}.sign(idx) < 0)));
    ki = kappa(IE, 1 + (s < 0));
    num = tau*(1 - kj)/dx + (al(ie)*delta_2/4 + N{iN}.sign(idx) .* uI(1, N{iN}.iE(idx))/2) .* (1 + kj);
    Mi(ie, idx) = num / (tau*(ki - 1)/dx + (al(ie)*delta_2/2 - s*uI(1, IE)/2) * (1 + ki));
end
rho = max(abs(eig(Mi)));
end


function kappa = update_kappa(kappa, alpha_opt, iN, root, uI, tau, dx, nx, N, delta_2)
%UPDATE_KAPPA  Transmission factors of the inner edges of node iN at their
%   OTHER end, determined by the Robin condition at node iN
%   (column 1: end x = 0, column 2: end x = 1).
for ie = 1:numel(N{iN}.iE)
    if N{iN}.ad(ie) > 0
        IE = N{iN}.iE(ie);  s = N{iN}.sign(ie);
        c = (-s*uI(1, IE) + alpha_opt(IE, iN)*delta_2) / 2;
        r1 = root(1, IE);  r2 = root(2, IE);
        if s < 0      % edge ends at iN (x = 1): factor at x = 0
            beta = (r2*(c - tau/dx) + c + tau/dx) / (r1*(c - tau/dx) + c + tau/dx);
            kappa(IE, 1) = (r2*r1^(nx-1) - beta*r1*r2^(nx-1)) / (r1^(nx-1) - beta*r2^(nx-1));
        else          % edge starts at iN (x = 0): factor at x = 1
            beta = (r2*(c + tau/dx) + c - tau/dx) / (r1*(c + tau/dx) + c - tau/dx);
            kappa(IE, 2) = (r2^(nx-1) - beta*r1^(nx-1)) / (r2^nx - beta*r1^nx);
        end
    end
end
end