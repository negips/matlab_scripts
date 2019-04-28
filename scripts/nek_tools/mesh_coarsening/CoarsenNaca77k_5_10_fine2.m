function ifc = CoarsenNaca77k_5_10_fine2(LX,LY,j,i,iflocked)

%    ARcut = 2.5;             % Coarsen if Aspect ratio is larger than this.
    skip_layers = 1;

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

    if i==28
        if j>=118 && j<=350
          ifc=1;
%        if xmid<=0.06
%          ifc=1;
        end  
    elseif i==30
%       Lower side          
        if j>=54 && j<=167
          ifc=1;
        end  
%       Upper side
        if j>=182 && j<=294
          ifc=1;
        end  
    elseif i==32
%       Lower side          
        if j>=58 && j<=84
          ifc=1;
        end  
%       Upper side
        if j>=151 && j<=174
          ifc=1;
        end  
    end  

%   Radially outgoing stuff      
%    if xmid<0.02
%      if i==start_layer+1
%        ifc=1;
%      end
%    end  

    
%   Skip first n layers      
    if i<=skip_layers
      ifc=0;
    end  

    [pts nels]=size(LX);  
%   End condition           
    if (iflocked(j) || j==1 || j==nels-1)
      ifc=0;
    end


end % function
%---------------------------------------------------------------------- 

