function [S, q, rho] = couplingm(example, g, nE, nN, N)
% Computes coupling conditions at network nodes.
% Solves a linear system based Half-moment approximation and kinetic coupling conditions
arguments
    example (1,1) double
    g (:, :) double
    nE (1,1) double
    nN (1,1) double
    N cell
end
switch example
    case 0
        del1 = 4/(sqrt(2*pi));
        alpha = (5*pi - 16) / (sqrt(2*pi) * (pi - 3));
        beta = (3*pi - 10) / (pi - 3);
        gamma = (3*pi - 8) / (sqrt(2*pi) * (pi - 3));
        lambda = sqrt(alpha + 2*gamma) / (nthroot(2*pi, 4) * alpha);

        % Block matrix A for global linear system
        b = zeros(nE * 4, 1);

        A = [0 0 (-1/sqrt(2*pi)) 0 1/2 0 (((gamma*beta - alpha)/(alpha*gamma) - lambda)/2) 0; (-1/(2*sqrt(2*pi))) 0 1/2 0  (-1/(2*sqrt(2*pi))) 0 -1/2 0; 1/2 0 (-2/sqrt(2*pi)) 0 0 0 -lambda*alpha/(2*gamma) 0; ...
            0 0 0 (-1/sqrt(2*pi)) 0 1/2 0 (((gamma*beta - alpha)/(alpha*gamma) - lambda)/2); 0 (-1/(2*sqrt(2*pi))) 0 1/2 0  (-1/(2*sqrt(2*pi))) 0 -1/2 ; 0 1/2 0 (-2/sqrt(2*pi)) 0 0 0 -lambda*alpha/(2*gamma);...
            0 0 1 1 0 0 0 0 ; 1 -1 del1 del1 0 0 0 0];
        b(1:3,1) = g(:,1);
        b(4:6,1) = g(:,2);

        % Solve linear system
        m = linsolve(A, b);

        % Return coupled variables
        S   = m(1:nE)';
        q   = m(nE+1 : 2*nE)';
        rho = m(nE*2+1 : nE*3)';

    otherwise

        del1 = 4/(3*sqrt(2*pi));

        % Block matrix A for global linear system
        A = zeros(nE * 3);
        b = zeros(nE * 3, 1);

        k = 0;
        j = 0;
        for iN = 1:nN
            for ie = 1:length(N{iN}.iE)
                IE = N{iN}.iE(ie);
                aD = N{iN}.ad(ie);

                % Edge equations (from kinetic equation)
                switch ie
                    case 1
                        A((iN-1)*3 + 2, IE) = 1;
                        A((iN-1)*3 + 1 : (iN-1)*3 + 2, IE+nE) = [1; del1];
                    case 2
                        A((iN-1)*3 + 2 : (iN-1)*3 + 3, IE) = [-1; 1];
                        A((iN-1)*3 + 1 : (iN-1)*3 + 3, IE+nE) = [1; -del1; del1];
                    case 3
                        A((iN-1)*3 + 3, IE) = -1;
                        A((iN-1)*3 + 1 : (iN-1)*3 + 3, IE+nE) = [1; 0; -del1];
                end

                if N{iN}.sign(ie) < 0
                    A((iN-1)*3+1:iN*3, IE+nE) = -A((iN-1)*3+1:iN*3, IE+nE);
                end

                % Handling boundary conditions (Half-moment-approximation)
                switch aD
                    case 0 % External boundary
                        k = k + 1;
                        b(nN*3 + (k-1)*2 + 1 : nN*3 + k*2) = g(:, IE);

                        if N{iN}.sign(ie) < 0
                            A(nN*3 + (k-1)*2 + 1 : nN*3 + k*2, IE) = [1/(2*sqrt(2*pi)); 1/2];
                            A(nN*3 + (k-1)*2 + 1 : nN*3 + k*2, nE+IE) = [1/2; 2/sqrt(2*pi)];
                            A(nN*3 + (k-1)*2 + 1 : nN*3 + k*2 - 1, nE*2+IE) = 1/(2*sqrt(2*pi));

                        else
                            A(nN*3 + (k-1)*2 + 1 : nN*3 + k*2, IE) = [-1/(2*sqrt(2*pi)); 1/2];
                            A(nN*3 + (k-1)*2 + 1 : nN*3 + k*2, nE+IE) = [1/2; -2/sqrt(2*pi)];
                            A(nN*3 + (k-1)*2 + 1, nE*2+IE) = -1/(2*sqrt(2*pi));

                        end
                    otherwise % Internal node (dummy equation for regularity reason)
                        if aD > iN
                            A(end - j, nE*2+IE) = 1;
                            j = j + 1;
                        end
                end
            end
        end
        % Solve linear system
        m = linsolve(A, b);

        % Return coupled variables
        S   = m(1:nE)';
        q   = m(nE+1 : 2*nE)';
        rho = m(nE*2+1 : nE*3)';
end
end