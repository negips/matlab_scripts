function [area] = areafcn(time,var,lag);

ae   = getalpha(time,lag);
area = abs(trapz(ae,var));

end
