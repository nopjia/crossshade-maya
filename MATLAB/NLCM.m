function [ x ] = NLCM( N, NC, M, MI)
% Solves a nonlinear constrained minimization problem
% Inputs:  N = number of cross section curves
%          NC = number of cross section intersections
%          M = matrix of cross section curve information
%          M = [t(i,j)_x t(i,j)_y x(i,j)_x x(i,j)_y ...
%                                      ]
%          MI = matrix of intersections, i.e if curves 1 and 3 intersect
%          then MI(1,3) == 1 else M(1,3) == 0
% Outputs: x = matrix of calculated normals and t_z 
%   Eric Lee

    % Constraints function
    function [c, ceq] = constraints(x)
        % Epsilon - value pulled from paper
        e = 0.1;
        
        % Nonlinear inequality constraints
        NIC = zeros(1, NC*4);
        c1 = 0;
        for k=1:N
            for l=k:N
                if k ~= l && MI(k,l) == 1
                    % Normals dot product
                    NIC(4*(c1) + 1) = (x(3*(k-1) + 1)*x(3*(l-1) + 1) + ...
                                       x(3*(k-1) + 2)*x(3*(l-1) + 2) + 1.0) - e;
                    
                    NIC(4*(c1) + 2) = -(x(3*(k-1) + 1)*x(3*(l-1) + 1) + ...
                                       x(3*(k-1) + 2)*x(3*(l-1) + 2) + 1.0) - e;
                    
                    % Tangents dot product
                    NIC(4*(c1) + 3) = (M(k,1)*M(l,1) + M(k,2)*M(l,2) + x(3*k)*x(3*l)) - e;
                    NIC(4*(c1) + 4) = -(M(k,1)*M(l,1) + M(k,2)*M(l,2) + x(3*k)*x(3*l)) - e;
                    c1 = c1 + 1;
                end
            end
        end
        c = NIC;
        c
        
        % Nonlinear equality constraints

        % Equation 5 from the paper
        % INSERT ACTUAL VALUES FOR T_X AND T_Y
        % Then equation 7 from the paper
        % AGAIN, INSERT ACTUAL VALUES FOR T_X AND T_Y      
        TN = zeros(1,N);
        for j=1:N
            TN(j) = M(j,1)*x(3*(j-1) + 1) + M(j,2)*x(3*(j-1) + 2) + x(3*j);
            %disp(3*(j));
        end
        ceq = TN;

    end
    
    % Function to be minimized
    function f = cmFunction(x)
        tempF = 0;
        c = 0;
        for k=1:N
            for l=k:N
                if k ~= l && MI(k,l) == 1
                    n0 = [x(3*(k-1) + 1) x(3*(k-1) + 2) 1];
                    n1 = [x(3*(l-1) + 1) x(3*(l-1) + 2) 1];
                    %t10 = [M(l,1), M(l,2), x(3*l)];
                    %t01 = [M(k,1), M(k,2), x(3*k)];
                    t10 = [M(2*c + 2,1), M(2*c + 2,2), x(3*l)];
                    t01 = [M(2*c + 1,1), M(2*c + 1,2), x(3*k)];

                    c = c + 1;
                    tempF = tempF + ((norm(cross(t10, n0),2)).^2 + ...
                        (norm(cross(t01, n1),2)).^2) + x(3*k).^2 + x(3*l).^2;     
                end
            end
        end
        
        f = tempF;
        %f
    end

    % Initial guess vector
    % n1, n2, t_z
    % BASED ON ORIENTATION CHOICE
    %x0 = zeros(1,N*3);
    %for i=1:(N*3) 
    x0 = zeros(1,N*2 + NC*2);
    for i=1:(N*2 + NC*2) 
        % 0 is initial guess for t_z
        %if mod(i,3) == 0 
         %   x0(i) = 0;
        % otherwise initial guess should be based on orientation choice
        %else
            x0(i) = 1;    
        %end
    end
    %x0 

    % Algorithm options
    options = optimset('Algorithm','interior-point');
    
    % params: function to minimize, initial guess, inequalities and bounds, 
    % constraints function
    [x,fval] = fmincon(@cmFunction,x0,[],[],[],[],[],[],@constraints,options)
    
    %Calculate offsets
    % Normals
    %n0 = [x(1) x(2) 1.0];
    %n1 = [x(3) x(4) 1.0];

    %solve(-9.4343563517316369*(n0(1) - n1(1)) + 1.3981787820170251*(n0(2) - n1(2)) + (0 - c), c)
end

