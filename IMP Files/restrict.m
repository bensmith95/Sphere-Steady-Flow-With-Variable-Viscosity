%% Restrict to Coarse grid
function r2h = restrict(rh)
    % Initialise
    N = 0.5*(size(rh)-3)+1; r2h = zeros(N+2);
    
    % Lagrange Interpolation
    for j = 3:N(1)
        for i = 3:N(2)
            r2h(j,i) = 0.0625*(rh(2*j-3,2*i-3)+rh(2*j-3,2*i-1)+rh(2*j-1,2*i-3)+rh(2*j-1,2*i-1)) + ...
                     0.125*(rh(2*j-2,2*i-3)+rh(2*j-2,2*i-1)+rh(2*j-3,2*i-2)+rh(2*j-1,2*i-2)) + ...
                     0.25*rh(2*j-2,2*i-2);
        end
    end
    % Injection
    r2h([2,N(1)+1],2:N(2)+1) = rh([2,end-1],2:2:end-1); r2h(2:N(1)+1,[2,N(2)+1]) = rh(2:2:end-1,[2,end-1]);
end