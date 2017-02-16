function val = barrierP3( x )
    % Barrier function for problem 3

    g1 = @(x) x(1) + x(2) - 3;
    g2 = @(x) -x(1) + 2.*x(2) - 4;

    if (x(1) + x(2) > 3) || (-x(1) + 2*x(2) > 4)
        val = 1e99;
    else
        val = -1/g1(x) -1/g2(x);
    end
end

