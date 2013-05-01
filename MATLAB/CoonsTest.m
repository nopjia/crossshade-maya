% Coons patch testing

clear all;
close all;
%{
X = [0 0.25 0.5 0.75 1
     0 0 0 0 1
     0 0 0 0 1
     0 0 0 0 1
     0 0.25 0.5 0.75 1];
 
 Y = [1 1 1 1 1
     0.75 0 0 0 0.75
     0.5 0 0 0 0.5
     0.25 0 0 0 0.25
     0 0 0 0 0];
 
 Z = [0 0.25 0.5 0.25 0
     0.25 0 0 0 0.25
     0.5 0 0 0 0.5
     0.25 0 0 0 0.25
     0 0.25 0.5 0.25 0];
 %}

X = [0 0.25 0.5 0.75 1
     0.2 0 0 0 1.2
     0.4 0 0 0 1.4
     0.2 0 0 0 1.2
     0 0.25 0.5 0.75 1];
 
 Y = [1 1.2 1.4 1.2 1
     0.75 0 0 0 0.75
     0.5 0 0 0 0.5
     0.25 0 0 0 0.25
     0 0.2 0.4 0.2 0];
 
 Z = [0 0.25 0.5 0.25 0
     0.25 0 0 0 0.25
     0.5 0 0 0 0.5
     0.25 0 0 0 0.25
     0 0.25 0.5 0.25 0];
 
 [m,n] = size(X);

% Iterate through matrix of points (x,y,z)
for u = 2:(m-1)
    for v = 2:(n-1)
        disp('Evaluating interior points');
        X(u,v) = coons(X,u,v,m,n);
        Y(u,v) = coons(Y,u,v,m,n);
        Z(u,v) = coons(Z,u,v,m,n);
    end
end

X
Y
Z

% Plot 3D surface
surf(X,Y,Z)