function [ x ] = twoCurveNLCM( A )
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
            % This test case works
             %(0.87091212321166278*0.97982164446094777 + 0.49143877914034567*-0.19987382281290467 + x(5)*x(6)) - e;
            %-(0.87091212321166278*0.97982164446094777 + 0.49143877914034567*-0.19987382281290467 + x(5)*x(6)) - e];
            (A(1,1)*A(2,1) + A(1,2)*A(2,2) + x(5)*x(6)) - e;
            -(A(1,1)*A(2,1) + A(1,2)*A(2,2) + x(5)*x(6)) - e];
        
        c


         %t_ij dot n_i
         ceq = [(A(1,1)*x(1) + A(1,2)*x(2) + x(5));
                (A(2,1)*x(3) + A(2,2)*x(4) + x(6))];
    end

    function f = cmFunction(x)
        n0 = [x(1) x(2) 1];
        n1 = [x(3) x(4) 1];
        
        t10 = [A(2,1), A(2,2), x(6)];
        t01 = [A(1,1), A(1,2), x(5)];

        f = ((norm(cross(t10, n0),2)).^2 + (norm(cross(t01, n1),2)).^2) + x(5).^2 + x(6).^2;
        %f
    end

x0 = [1,1, 1,1, 1,1]; 

options = optimset('Algorithm','interior-point');

[x,fval] = fmincon(@cmFunction,x0,[],[],[],[],[],[],@constraints,options)

end

