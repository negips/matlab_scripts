function [lambda]=SimsonStep(x,xmin,xmax,xrise,xfall)

   lambda_max = 1;   
   lambda = lambda_max*(Step((x-xmin)/xrise) - Step((x-xmax)/xfall + 1));

   return

function y = Step(x)

  l=length(x);
  ind1=find(x<=0);
  ind2=find(x>=1);
  ind3=1:l;
  ind3([ind1; ind2])=[];
      
  y=x;
  y(ind1)=0;
  y(ind2)=1;
  x2=x(ind3)
  if ~(isempty(x2))    
    y(ind3)=1./( 1 + exp( (1/(x2-1) + 1/x2) ) );
  end

%  if (max(y)>1)
%    dbstop in SimsonStep at 27
%  end         
%  ;    

  return     
      
