function pts = georefine(pt0,pt1,Npt,pratio)
            
            if pratio == 1.0
                pratio = 1.0+1.0e-10;
            end
            
            pts = zeros(1,Npt);
            pts(1) = pt0;
            
            fc  = (pt1-pt0)*(1.0-pratio)/(1.0-pratio^(Npt-1));
            for p = 2:Npt
                pts(p) = pts(p-1)+fc*pratio^(p-2);
            end
            
end
