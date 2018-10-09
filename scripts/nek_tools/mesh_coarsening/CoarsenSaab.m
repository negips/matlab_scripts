function ifc = CoarsenSaab(LX,LY,j,i,iflocked)

    ARcut = 2.5;             % Coarsen if Aspect ratio is larger than this.
    start_layer = 13;

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

    xmid = mean(LX(:,j));
    ymid = mean(LY(:,j));
    rad = sqrt(xmid^2 + ymid^2);  

    ifc=0;

    if ymid>0.00
%     Upper Side          
      if xmid>0.04
        if i==start_layer+1
          ifc=1;
        elseif i==start_layer+3
          ifc=1;
        elseif i==start_layer+5
          ifc=1;
        else    
          ifc=0;
        end
      end
      
    else 
%     Lower Side          

    end  % ymid>0

%   Radially outgoing stuff      
%    if xmid<0.02
%      if i==start_layer+1
%        ifc=1;
%      end
%    end  

    
%   Skip first n layers      
    if i<start_layer
      ifc=0;
    end  

    [pts nels]=size(LX);  
%   End condition           
    if (iflocked(j) || j==1 || j>=nels-1)
      ifc=0;
    end


end % function
%---------------------------------------------------------------------- 

