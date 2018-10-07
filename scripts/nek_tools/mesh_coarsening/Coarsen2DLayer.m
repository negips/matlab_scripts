function [LX,LY,ifc,NewX,NewY,NewCEl,NewCoF,NewET] = Coarsen2DLayer(LE,LX,LY,ifc,NewX,NewY,NewCEl,NewCoF,NewET,i,fig2,ifplot)

  l1 =length(LX);
  cmap = jet(l1);
  skipnext=0;
  skipnext2=0;
  skipnext3=0;
  ifcontinue=0; 
  for j=1:l1

%   Skipping coarsening for the first element
    if j==1
      ifcontinue=1;
    end  
      
%   Skip subsequent elements after coarsening  
    if skipnext
      skipnext=0;
      ifcontinue=1;
    elseif skipnext2
      skipnext2=0;
      ifcontinue=1;
    elseif skipnext3
      skipnext3=0;
      ifcontinue=1;
    end

%   If not coarsen skip to next iteration
    if (~ifc(j))
      ifcontinue=1;
    end  

    if (ifcontinue)
      if (ifplot)
        figure(fig2)
        xt = LX(:,j);
        yt = LY(:,j);
        figure(2)
        fill(xt,yt,cmap(j,:)); hold on
      end  
      
      ifcontinue=0;
      continue
    end   

%   Coarsen Element j
    LX(4,j)=NewX{i+1}(4,j);
    LY(4,j)=NewY{i+1}(4,j);
%   Connecting Element no changes  
    NewCEl{i}(4,j)=NewCEl{i}(4,j-1);  
%   onFace of connecting element changes
    NewCoF{i}(4,j)=3;
%   Element type
    NewET{i}{j}='e4';

    if j<l1
%     Coarsen Element j+1
      LX(1,j+1)=LX(4,j);
      LY(1,j+1)=LY(4,j);
%     Connecting Element no changes  
      NewCEl{i}(4,j+1)=NewCEl{i}(4,j+2); 
%     onFace of connecting element changes
      NewCoF{i}(4,j+1)=1;
%     Element type      
      NewET{i}{j+1}='e4';
    end  

    if (ifplot)
      figure(fig2)
      xt = LX(:,j);
      yt = LY(:,j);
      figure(2)
      fill(xt,yt,cmap(j,:)); hold on
    end

%   Skip next few elements
    skipnext=1;
    skipnext2=1;
    skipnext3=1;
    if (j+1<=l1)
      ifc(j+1)=0;
    end  
    if (j+2<=l1)
      ifc(j+2)=0;
    end  
    if (j+3<=l1)
      ifc(j+3)=0;
    end 

  end % end j

  NewX{i}=LX;
  NewY{i}=LY; 

end   % end function

%---------------------------------------------------------------------- 

