function ifc = CoarsenCriteria(LX,LY,j,i,iflocked)

    ARcut = 2.5;             % Coarsen if Aspect ratio is larger than this.
    start_layer = 5;

%   Find aspect ratio of elements 

%   Length of sides 1 and 3   % Facing 'O'
    l1o = sqrt( (LX(1,j)-LX(2,j))^2 + (LY(1,j) - LY(2,j))^2);
    l2o = sqrt( (LX(3,j)-LX(4,j))^2 + (LY(3,j) - LY(4,j))^2);
    dlo = mean([l1o l2o]);
     
%   Length of sides 2 and 4   % Facing 'V'
    l1v = sqrt( (LX(2,j)-LX(3,j))^2 + (LY(2,j) - LY(3,j))^2);
    l2v = sqrt( (LX(4,j)-LX(1,j))^2 + (LY(4,j) - LY(1,j))^2);
    dlv = mean([l1v l2v]);

    l_ar = dlo/dlv;
    lmax = max([l1o l2o l1v l2v]);

    ifc=0;
    if l_ar>ARcut
      ifc=1; 
    end

    xmid = mean(LX(:,j));
    ymid = mean(LY(:,j));
    rad = sqrt(xmid^2 + ymid^2);  

    if (xmid<0.0 && rad>0.1 )
      ifc=0;
    end

    if xmid>1 && l_ar>1.25
      ifc=1;
    end

%   Maximum length  
    if dlv>0.1
      ifc=0;
    end

%   For the radially emerging elements I refine by number of layers
    if i==start_layer
      theta=atan(ymid/(xmid-0.25))*180/pi;
      if xmid<0.25 && abs(theta)<15
        ifc=1;
      end
    elseif i<=start_layer+4
      theta=atan(ymid/(xmid-0.25))*180/pi;
      if xmid<0.25 && abs(theta)>15 && abs(theta)<75
        ifc=1;
      end
    else
      if xmid<0.25 
        ifc=0;
      end
    end
    
%   Skip first n layers      
    if i<start_layer
      ifc=0;
    end  

    [pts nels]=size(LX);  
%   End condition           
    if (iflocked(j) || j==1 || j>=nels-1)
      ifc=0;
    end
    
%   If this layer has anything locked, lock the whole layer
    if max(iflocked)
      ifc=0;
    end  


end % function
%---------------------------------------------------------------------- 

