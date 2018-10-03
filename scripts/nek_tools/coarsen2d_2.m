% First attempts at coarsening the mesh

clear
clc
close all

load saab_wing2d.mat

skiplayers = 10;         % No of layers to skip when coarsening
lmax = 0.1;              % Maximum length of a side. If the length is larger. Don't coarsen
ARcut = 2.0;             % Coarsen if Aspect ratio is larger than this.


for i=1:nlayers
  j=nlayers-i+1;
  OldE{i}=LayerE{j};
  OldX{i}=LayerX{j};
  OldY{i}=LayerY{j};
  OldBC{i}=LayerBC{j};

  NewE{i}=OldE{i};
  NewX{i}=OldX{i};
  NewY{i}=OldY{i};
  NewBC{i}=OldBC{i};
end


% Test coarsening algorithm 1
% Coarsen Layer by Layer
% Define aspect ratio as 'O' face lengths to 'V' face lengths
layer_start = skiplayers+1;
for i=1:layer_start  %nlayers

  if (i<layer_start)
    continue
  end  

  LE=NewE{i};
  LX=NewX{i};
  LY=NewY{i};

  l1 =length(LE);
  cmap = jet(l1); 
  for j=1:l1

%   Find aspect ratio of elements 

%   Length of sides 1 and 3   % Facing 'O'
    l1 = sqrt( (LX(1,j)-LX(2,j))^2 + (LY(1,j) - LY(2,j))^2);
    l2 = sqrt( (LX(3,j)-LX(4,j))^2 + (LY(3,j) - LY(4,j))^2);
    dlo = mean([l1 l2]);
     
%   Length of sides 2 and 4   % Facing 'V'
    l1 = sqrt( (LX(2,j)-LX(3,j))^2 + (LY(2,j) - LY(3,j))^2);
    l2 = sqrt( (LX(4,j)-LX(1,j))^2 + (LY(4,j) - LY(1,j))^2);
    dlv = mean([l1 l2]);

    l_ar(j) = dlo/dlv;

    xt = LX(:,j);
    yt = LY(:,j);

    figure(1)
    fill(xt,yt,cmap(j,:)); hold on

    ifc(j)=0;
    if l_ar(j)>ARcut
      ifc(j)=1; 
    end  
  end 

  c_ind = find(ifc);
  nc = length(c_ind);
  disp(['Coarsen ' num2str(nc) ' Elements in layer ' num2str(i)])
  for j=1:nc
    k=c_ind(j);
    xmid=mean(LX(:,k));
    ymid=mean(LY(:,k));
    zmid=2;
    figure(1);
    plot3(xmid,ymid,zmid, '. ', 'MarkerSize', 12);
  end

% Coarsen layer in consecutive pairs
%  coarsen_layer
  if nc>0
    fig2=figure(2);    
    [LX,LY,ifc,NewX,NewY] = CoarsenLayer(LE,LX,LY,ifc,NewX,NewY,i,fig2);
%   Modify all subsequent Layers
    if3skip=1;  
    [NewE, NewX, NewY, NewBC]=CreateNewLayers(NewE,NewX,NewY,NewBC,i,ifc,if3skip,fig2);
  end

end

%---------------------------------------------------------------------- 

function [LX,LY,ifc,NewX,NewY] = CoarsenLayer(LE,LX,LY,ifc,NewX,NewY,i,fig2)

  l1 =length(LE);
  cmap = jet(l1);
  skipnext=0;
  skipnext2=0;
  skipnext3=0;
  incontinue=0; 
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
      figure(fig2)
      xt = LX(:,j);
      yt = LY(:,j);
      figure(2)
      fill(xt,yt,cmap(j,:)); hold on
      
      ifcontinue=0;
      continue
    end   

    LX(4,j)=NewX{i+1}(4,j);
    LY(4,j)=NewY{i+1}(4,j);

%   Next Element
    LX(1,j+1)=LX(4,j);
    LY(1,j+1)=LY(4,j);

    figure(fig2)
    xt = LX(:,j);
    yt = LY(:,j);
    figure(2)
    fill(xt,yt,cmap(j,:)); hold on

%   Skip next element    
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

function [NewE, NewX, NewY, NewBC]=CreateNewLayers(NewE,NewX,NewY,NewBC,cr_layer,ifc,if3skip,fig2)
 
  nlayers=length(NewX);
  il=0;
  for i=cr_layer+1:nlayers
    il=il+1;  
    LE=NewE{i};
    LX=NewX{i};
    LY=NewY{i};
    
    LX2=LX;
    LY2=LY;
    l1=length(LE);


    ifcontinue=0;
    k=0;          % Index for LX2, LY2
    for j=1:l1

      k=k+1;

      if (~ifc(j))
        continue
      end  

      if (il==1)  % First layer. Needs slightly different treatment

%       Enlarge K-1th element            
        LX2(1,k-1)=LX(1,j-1); 
        LX2(2,k-1)=LX(2,j-1); 
        LX2(3,k-1)=LX(3,j-1); 
        LX2(4,k-1)=LX(4,j); 
   
        LY2(1,k-1)=LY(1,j-1); 
        LY2(2,k-1)=LY(2,j-1); 
        LY2(3,k-1)=LY(3,j-1); 
        LY2(4,k-1)=LY(4,j);

%       Enlarg K+2th element        
        LX2(1,k+2)=LX(1,j+1); 
        LX2(2,k+2)=LX(3,j+1); 
        LX2(3,k+2)=LX(3,j+2); 
        LX2(4,k+2)=LX(4,j+2); 
   
        LY2(1,k+2)=LY(1,j+1); 
        LY2(2,k+2)=LY(3,j+1); 
        LY2(3,k+2)=LY(3,j+2); 
        LY2(4,k+2)=LY(4,j+2);

%       Delete K and K+1th element
        LX2(:,[k k+1]) = [];
        LY2(:,[k k+1]) = [];
      else 

%       Enlarge K-1th element            
        LX2(1,k-1)=LX(1,j-1); 
        LX2(2,k-1)=LX(2,j-1); 
        LX2(3,k-1)=LX(3,j); 
        LX2(4,k-1)=LX(4,j); 
   
        LY2(1,k-1)=LY(1,j-1); 
        LY2(2,k-1)=LY(2,j-1); 
        LY2(3,k-1)=LY(3,j); 
        LY2(4,k-1)=LY(4,j);

%       Enlarg K+2th element        
        LX2(1,k+2)=LX(1,j+1); 
        LX2(2,k+2)=LX(2,j+1); 
        LX2(3,k+2)=LX(3,j+2); 
        LX2(4,k+2)=LX(4,j+2); 
   
        LY2(1,k+2)=LY(1,j+1); 
        LY2(2,k+2)=LY(2,j+1); 
        LY2(3,k+2)=LY(3,j+2); 
        LY2(4,k+2)=LY(4,j+2);

%       Delete K and K+1th element
        LX2(:,[k k+1]) = [];
        LY2(:,[k k+1]) = [];

      end

      k=k-2;

    end 

    NewX{i}=LX2;
    NewY{i}=LY2;

    l1=length(LX2);
    cmap=jet(l1);
    for j=1:l1
      xt=LX2(:,j);
      yt=LY2(:,j);
      figure(fig2);
      fill(xt,yt,cmap(j,:)); hold on
    end
  end  

end   % end function

%---------------------------------------------------------------------- 





