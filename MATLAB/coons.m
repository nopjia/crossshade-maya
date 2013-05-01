function [ result ] = coons( M, u, v, degreeU, degreeV)
% Evaluates Coons patch surface parameterization function 
%   

% Bottom to top interpolation  
U = double(u);
mU = double(degreeU);
  
uOvermU = U/mU;
value1 = (1.0 - uOvermU) * M(1,v) + (uOvermU * M(mU,v));

% Left to right interpolation
V = double(v);
mV = double(degreeV);
vOvermV = V/mV;
value2 = (1.0 - vOvermV) * M(u,1) + (vOvermV * M(u,mV));

% Bilinear interpolation
a = M(1,1); % lower left corner
b = M(mU,1); % upper left corner
c = M(1,mV); % lower right corner
d = M(mU,mV); % upper right corner
value3 = (a+(b-a)*uOvermU) + ((c+(d-c)*uOvermU)-(a+(b-a)*uOvermU))*vOvermV;


result = value1 + value2 - value3;
%{
value = (2 - u)*S(1,v) + u*S(2,v) + ...
    (2 - v)*S(u,1) + v*S(u,2) - ...
    [2 - u u] * M *M2;
%}

end

