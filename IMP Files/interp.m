%% Interpolate to Fine grid
function F = interp(C)
    % Initialise
    N = size(C); F = zeros(2*(N-3)+1); M = size(F)+2; F = zeros(M);
    
    % Inject
    F(2:2:M(1)-1,2:2:M(2)-1) = C(2:N(1)-1,2:N(2)-1);
    
    % Linear Interpolation/Average
    F(2:2:M(1)-1,3:2:M(2)-2) = 0.5*(C(2:N(1)-1,2:N(2)-2)+C(2:N(1)-1,3:N(2)-1));
    F(3:2:M(1)-2,2:2:M(2)-1) = 0.5*(C(2:N(1)-2,2:N(2)-1)+C(3:N(1)-1,2:N(2)-1));  
    F(3:2:M(1)-2,3:2:M(2)-2) = 0.25*(C(2:N(1)-2,2:N(2)-2)+C(2:N(1)-2,3:N(2)-1) + ...
                                     C(3:N(1)-1,2:N(2)-2)+C(3:N(1)-1,3:N(2)-1));
end