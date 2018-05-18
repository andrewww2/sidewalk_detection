% Normalize matrix
function[A] = NormalizeMatrix(A)
    A = A/max(A(:));
end