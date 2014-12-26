
function dx = f(t,x);
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = -sin(x(1)); %y'' = - sin(y)
