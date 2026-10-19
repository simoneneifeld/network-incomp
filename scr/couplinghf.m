function [S, q, rho] = couplinghf(example, g, nE, nN, N)
% Computes coupling conditions at network nodes.
% Solves a linear system based on half-moment approximation and kinetic coupling conditions
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
        del1 = 0.528767;  %0.5298; %0.528767
        alpha = (5*pi - 16) / (sqrt(2*pi) * (pi - 3));
        beta = (3*pi - 10) / (pi - 3);
        gamma = (3*pi - 8) / (sqrt(2*pi) * (pi - 3));
        lambda = sqrt(alpha + 2*gamma) / (nthroot(2*pi, 4) * alpha);

        A = zeros(nE * 4);
        b = zeros(nE * 4, 1);
        k = 0;
        j = 0;

        % Block matrix A for global linear system
        for iN = 1:nN
            for ie = 1:length(N{iN}.iE)
                IE = N{iN}.iE(ie);
                aD = N{iN}.ad(ie);

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

                switch aD
                    case 0
                        k = k + 1;
                        b(nN*3 + (k-1)*3 + 1 : nN*3 + k*3) = g(:, IE);
                        if N{iN}.sign(ie) < 0
                            A(nN*3 + (k-1)*3 + 2 : nN*3 + k*3, IE) = [1/(2*sqrt(2*pi)); 1/2];
                            A(nN*3 + (k-1)*3 + 1 : nN*3 + k*3, nE+IE) = [1/sqrt(2*pi); 1/2; 2/sqrt(2*pi)];
                            A(nN*3 + (k-1)*3 + 1 : nN*3 + k*3 - 1, nE*2+IE) = [1/2; 1/(2*sqrt(2*pi))];
                            A(nN*3 + (k-1)*3 + 1 : nN*3 + k*3, nE*3+IE) = [((gamma*beta - alpha)/(alpha*gamma) - lambda)/2; 1/2; -lambda*alpha/(2*gamma)];
                        else
                            A(nN*3 + (k-1)*3 + 2 : nN*3 + k*3, IE) = [-1/(2*sqrt(2*pi)); 1/2];
                            A(nN*3 + (k-1)*3 + 1 : nN*3 + k*3, nE+IE) = [-1/sqrt(2*pi); 1/2; -2/sqrt(2*pi)];
                            A(nN*3 + (k-1)*3 + 1 : nN*3 + k*3 - 1, nE*2+IE) = [1/2; -1/(2*sqrt(2*pi))];
                            A(nN*3 + (k-1)*3 + 1 : nN*3 + k*3, nE*3+IE) = [((gamma*beta - alpha)/(alpha*gamma) - lambda)/2; -1/2; -lambda*alpha/(2*gamma)];
                        end
                    otherwise
                        if aD > iN
                            A(end - j*2 - 1, nE*2+IE) = 1;
                            A(end - j*2, nE*3+IE) = 1;
                            j = j + 1;
                        end
                end
            end
        end
        m = A \ b;
        
        % Return coupled variables
        S   = m(1:nE)';
        q   = m(nE+1 : 2*nE)';
        rho = m(nE*2+1 : nE*3)';
end
end