% Testing out symbolic function in matlab

clear
clc

Tau = sym('Tau%d%d', [3 3]);
tau = sym('tau%d%d', [3 3]);

syms eta P p

Rot = [cos(eta) -sin(eta) 0;
       sin(eta) cos(eta)  0;
       0          0       1];

pos = sym('R%d', [3, 1]);
n   = sym('N%d', [3,1]);

Tau_force = transpose(Tau)*n;

C0 = cross(pos,Tau_force);

cm_P0 = cross(pos,P*n);

new_pos = Rot*pos;
new_n = Rot*n;
Tt = Tau + tau;
Tt_force = transpose(Tt)*new_n;

Ct = cross(new_pos,Tt_force);
cm_p = cross(new_pos,(P+p)*new_n);

Cex = expand(Ct-C0);

% Linearize deformation
C2 = subs(Cex,sin(eta),eta);
C3 = subs(C2,cos(eta),1);
C4 = subs(C3,n(3),0);
% Neglect higher order terms
C5 = subs(C4,eta^2,0);
for i=1:3
  for j=1:3
    C5 = subs(C5,tau(i,j)*eta,0);
  end
end  

coeff_eta = collect(C5(3),eta)


% For pressure
cm_pex = expand(cm_p-cm_P0);

% Linearize deformation
p2 = subs(cm_pex,sin(eta),eta);
p3 = subs(p2,cos(eta),1);
p4 = subs(p3,n(3),0);
% Neglect higher order terms
p5 = subs(p4,eta^2,0);
for i=1:3
  for j=1:3
    p5 = subs(p5,p*eta,0);
  end
end  

coeff_p = collect(p5(3),p)

% Testing using split sin and cos matrices

Cs  = [1 0 0;
       0 1 0;
       0 0 1];

Sn  = [0 -eta 0;
       eta 0  0;
       0   0  0];


% After structural linearization and dropping nonlinear terms

f_mean = cross(pos,transpose(Tau)*n);
f_eta1 = cross(pos,transpose(Tau)*(Sn*n));
f_eta2 = cross(Sn*pos,transpose(Tau)*n);
f_u    = cross(pos,transpose(tau)*n); 

f_perturb = f_eta1 + f_eta2 + f_u;
f_perturb2 = subs(f_perturb,n(3),0);
f_perturb3 = expand(f_perturb2);
coeff_eta2 = collect(f_perturb3(3),eta)





