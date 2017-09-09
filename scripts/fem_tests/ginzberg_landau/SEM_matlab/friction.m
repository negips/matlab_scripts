%friction.m	Friction coefficient
% mu = friction(u,f)
%   mu 	friction coefficient
%   u	slip
%   f	friction data
function mu = friction(u,f)
% Linear slip weakening friction law
  mu = max(f.MUs - f.W .* u ,f.MUd);

%-- Chambon's non linear law --
 % p = 0.4;
 % u = u ./(p*f.Dc);
 % mu = f.MUd +(f.MUs-f.MUd) ./ (1+u).^p ; 
% With this definition:
% initial slope = (MUs-MUd)/Dc
% if p=0.4 : 90% dynamic stress drop at slip = 126*Dc
%            99%  ...                        = 4e4*Dc
% so if Dc=100e-6 m: 90% drop at slip = 1.2 cm
%                    99% ...          = 4 m
