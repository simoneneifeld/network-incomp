function [nE, nN, N, pB, uB, rhoB] = setup_network_topology(example)
% Initializes network topology and macro boundary conditions
% based on the provided 'example'.

if (example==0)
    nE=2;
    nN=1;

    N{1}.iE = [1 2];
    N{1}.sign = [1 1];
    N{1}.ad = [0 0];

    pB = [1 2];
    uB = [0 0];
    rhoB = [1 1];

elseif (example == 1)
    nE = 3;
    nN = 1;

    N{1}.iE = [1 2 3];
    N{1}.sign = [1 1 1]; % 1 outgoing; -1 ingoing
    N{1}.ad = [0 0 0];
    
    %pB=[1 2 3];
    %uB = [0 0 0];
    %rhoB = [0 0 0];
    pB = [1 1 2];   % all values required
    uB = [0 0 0];
    rhoB = [1 1 1];

elseif (example == 2)
    nE = 7;
    nN = 4;

    N{1}.iE = [1 2 3];    N{1}.sign = [-1 1 1];    N{1}.ad = [0 2 3];
    N{2}.iE = [2 4 5];    N{2}.sign = [-1 1 1];    N{2}.ad = [1 3 4];
    N{3}.iE = [3 4 6];    N{3}.sign = [-1 -1 1];   N{3}.ad = [1 2 4];
    N{4}.iE = [5 6 7];    N{4}.sign = [-1 -1 1];   N{4}.ad = [2 3 0];

    pB = [1 0 0 0 0 0 2]; % only 1st and last value required
    uB = [0 0 0 0 0 0 0];
    rhoB = [1 0 0 0 0 0 1];

elseif (example == 3)
    nE = 5;
    nN = 2;

    N{1}.iE = [1 2 3];    N{1}.sign = [-1 -1 1];   N{1}.ad = [0 2 0];
    N{2}.iE = [2 4 5];    N{2}.sign = [1 -1 1];    N{2}.ad = [1 0 0];

    pB = [1 0 2 1 2]; % all except 2
    uB = [0 0 0 0 0];
    rhoB = [1 0 1 1 1];

elseif (example == 4)
    nE = 10;
    nN = 6;

    N{1}.iE = [1 2 3];    N{1}.sign = [-1 1 1];    N{1}.ad = [0 2 3];
    N{2}.iE = [2 4 5];    N{2}.sign = [-1 1 1];    N{2}.ad = [1 3 4];
    N{3}.iE = [3 4 6];    N{3}.sign = [-1 -1 1];   N{3}.ad = [1 2 5];
    N{4}.iE = [5 7 8];    N{4}.sign = [-1 -1 1];   N{4}.ad = [2 5 6];
    N{5}.iE = [6 7 9];    N{5}.sign = [-1 1 1];    N{5}.ad = [3 4 6];
    N{6}.iE = [8 9 10];   N{6}.sign = [-1 -1 1];   N{6}.ad = [4 5 0];

    pB = [1 0 0 0 0 0 0 0 0 2];
    uB = [0 0 0 0 0 0 0 0 0 0];
    rhoB = [1 0 0 0 0 0 0 0 0 1];
elseif (example == 5)
    nE = 8;
    nN = 4;

    N{1}.iE = [1 2 3];    N{1}.sign = [-1 1 1];    N{1}.ad = [0 2 3];
    N{2}.iE = [2 4 7];    N{2}.sign = [-1 -1 1];   N{2}.ad = [1 0 4];
    N{3}.iE = [3 5 6];    N{3}.sign = [-1 1 1];   N{3}.ad = [1 0 4];
    N{4}.iE = [6 7 8];    N{4}.sign = [-1 -1 1];   N{4}.ad = [3 2 0];

    pB = [1 0 0 1 2 0 0 2]; % only 1st and last value required
    uB = [0 0 0 0 0 0 0 0];
    rhoB = [1 0 0 1 1 0 0 1];
else
    error('Unknown example. Choose 0, 1, 2, 3, 4 or 5.');
end
end