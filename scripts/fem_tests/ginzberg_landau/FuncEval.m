function [val] = FuncEval(coeffs,x,ifderiv)
%    Evaluate Polynomial function
%    Or derivative at x

     n = length(coeffs);
     N = n-1;
     x_p_vec=zeros(n,1);

     if (ifderiv)                       % first entry is always zero.
          for p = 1:N
               x_p_vec(p+1) = x^(p-1);
          end
          f_p_vec = coeffs.*x_p_vec;
     else
          for p = 0:N
               x_p_vec(p+1) = x^p;
          end
          f_p_vec = coeffs.*x_p_vec;
     end
     val = sum(f_p_vec);
     return

