% spectral adjacency embedding
function [ASE] = adjSpecEmbed(A, d)
    A = diagAug(A);
    [V, D] = eigs(A, d);
    ASE = V * sqrt(abs(D));
end