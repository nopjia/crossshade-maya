function [ x ] = threeCurveNLCM( A )
% NLCM solver for two curve case
%   Eric Lee

    function [c, ceq] = constraints(x)
        % Epsilon - value pulled from paper
        e = 0.1;

        % Nonlinear inequality constraints
        %c = -x(1)*x(2) - 10; 
        c = [(x(1)*x(3) + x(2)*x(4) + 1.0) - e;
            -(x(1)*x(3) + x(2)*x(4) + 1.0) - e;
            %t_ij dot t_ji         
            (A(1,1)*A(2,1) + A(1,2)*A(2,2) + x(7)*x(8)) - e;
            -(A(1,1)*A(2,1) + A(1,2)*A(2,2) + x(7)*x(8)) - e;
            
            (x(3)*x(5) + x(4)*x(6) + 1.0) - e;
            -(x(3)*x(5) + x(4)*x(6) + 1.0) - e;
            (A(3,1)*A(4,1) + A(3,2)*A(4,2) + x(9)*x(10)) - e;
            -(A(3,1)*A(4,1) + A(3,2)*A(4,2) + x(9)*x(10)) - e];
        
        c


         %t_ij dot n_i
         ceq = [(A(1,1)*x(1) + A(1,2)*x(2) + x(7));
                (A(2,1)*x(3) + A(2,2)*x(4) + x(8));
                (A(3,1)*x(3) + A(3,2)*x(4) + x(9));
                (A(4,1)*x(5) + A(4,2)*x(6) + x(10))];
    end

    function f = cmFunction(x)
        n0 = [x(1) x(2) 1];
        n1 = [x(3) x(4) 1];
        
        t10 = [A(2,1), A(2,2), x(8)];
        t01 = [A(1,1), A(1,2), x(7)];

        s = ((norm(cross(t10, n0),2)).^2 + (norm(cross(t01, n1),2)).^2) + x(7).^2 + x(8).^2;
        
        n2 = [x(5) x(6) 1];
        
        t21 = [A(4,1) A(4,2) x(10)];
        t12 = [A(3,1) A(3,2) x(9)];
        
        f = s + ((norm(cross(t21, n1),2)).^2 + ...
            (norm(cross(t12, n2),2)).^2) + x(9).^2 + x(10).^2;
        
        %f
    end

x0 = [-1,-1, -1,-1, -1,-1  -1,-1 -1,-1]; 

options = optimset('Algorithm','interior-point');

[x,fval] = fmincon(@cmFunction,x0,[],[],[],[],[],[],@constraints,options)

end

