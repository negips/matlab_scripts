function [x pratio] = findgeoP(xs,xe,pratioinit,N,dx)
            
            %define anonymous function
            f = @(pratio)findP(xs,xe,pratio,N,dx);
            
            %find optimal pratio2
            options = optimset('GradObj','off','Display','off','LargeScale','off');
            [pratio,fval] = fminunc(f,pratioinit,options);
            
            
            x =  georefine(xs,xe,N,pratio);
            
            
            function f = findP(x0,x1,pratio,N,dxgiven)
                
                xd = georefine(x0,x1,N,pratio);
                
                dx = xd(2)-xd(1);
                
                f = abs(dxgiven - dx);
            end
            
        end