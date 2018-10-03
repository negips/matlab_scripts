% First attempts at coarsening the mesh

clear
clc
close all

load saab_wing2d.mat

skiplayers = 10;         % No of layers to skip when coarsening

disp(['Total Number of Elements: ', num2str(rea.mesh.nelg)])

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
new_nelg = 0;
fig2=figure(2);
for i=1:nlayers

  LE=NewE{i};
  LX=NewX{i};
  LY=NewY{i};    

  if (i<layer_start)
    nels_layer = length(LX);    
    cmap = jet(nels_layer);
    for j=1:nels_layer 
      xt = LX(:,j);
      yt = LY(:,j);

      figure(2)
      fill(xt,yt,cmap(j,:)); hold on
    end 
    
    new_nelg = new_nelg+nels_layer; 
    continue
  end  

  l1 =length(LX);
  if i==layer_start
    iflocked=zeros(l1,1);
  end  
  cmap = jet(l1);
  l2 = length(OldX{i});
  ifc = zeros(l2,1); 
  for j=1:l1
      
    ifc(j) = CoarsenCriteria(LX,LY,j,i,iflocked);
        
    xt = LX(:,j);
    yt = LY(:,j);

%    figure(1)
%    fill(xt,yt,cmap(j,:)); hold on
  end

  if i==nlayers
    ifc(:)=0;
  end  

  c_ind = find(ifc);
  nc = length(c_ind);
  disp(['Coarsen ' num2str(nc) ' Elements in layer ' num2str(i)])
%  for j=1:nc
%    k=c_ind(j);
%    xmid=mean(LX(:,k));
%    ymid=mean(LY(:,k));
%    zmid=2;
%    figure(1);
%    plot3(xmid,ymid,zmid, '. ', 'MarkerSize', 12);
%  end

%   Coarsen layer in consecutive pairs
    figure(2);    
    [LX,LY,ifc,NewX,NewY] = CoarsenLayer(LE,LX,LY,ifc,NewX,NewY,i,fig2);
%   Modify all subsequent Layers
    if3skip=1;  
    [NewE, NewX, NewY, NewBC, iflocked]=CreateNewLayers(NewE,NewX,NewY,NewBC,i,ifc,if3skip,fig2);
    nels_layer = length(NewX{i});
    new_nelg = new_nelg+nels_layer;
end
disp(['Total Number of Elements: ', num2str(new_nelg)])

%---------------------------------------------------------------------- 

function ifc = CoarsenCriteria(LX,LY,j,i,iflocked)

    ARcut = 2;             % Coarsen if Aspect ratio is larger than this.

%   Find aspect ratio of elements 

%   Length of sides 1 and 3   % Facing 'O'
    l1 = sqrt( (LX(1,j)-LX(2,j))^2 + (LY(1,j) - LY(2,j))^2);
    l2 = sqrt( (LX(3,j)-LX(4,j))^2 + (LY(3,j) - LY(4,j))^2);
    dlo = mean([l1 l2]);
     
%   Length of sides 2 and 4   % Facing 'V'
    l1 = sqrt( (LX(2,j)-LX(3,j))^2 + (LY(2,j) - LY(3,j))^2);
    l2 = sqrt( (LX(4,j)-LX(1,j))^2 + (LY(4,j) - LY(1,j))^2);
    dlv = mean([l1 l2]);

    l_ar = dlo/dlv;

    ifc=0;
    if l_ar>ARcut
      ifc=1; 
    end

    xmid = mean(LX(:,j));
    if (xmid<0.25)
      ifc=0;
    end

    if xmid>1 && l_ar>1.25
      ifc=1;
    end

    [pts nels]=size(LX);  
%   End condition           
    if (iflocked(j) || j==1 || j>=nels-1)
      ifc=0;
    end
     


end % function
%---------------------------------------------------------------------- 

function [LX,LY,ifc,NewX,NewY] = CoarsenLayer(LE,LX,LY,ifc,NewX,NewY,i,fig2)

  l1 =length(LX);
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

function [NewE, NewX, NewY, NewBC, iflocked]=CreateNewLayers(NewE,NewX,NewY,NewBC,cr_layer,ifc,if3skip,fig2)

  iflocked = [];

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

    if il==1
      iflocked=zeros(l1,1);
    end  

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

        iflocked(k-1)=1;
 
%       Enlarg K+2th element        
        LX2(1,k+2)=LX(1,j+1); 
        LX2(2,k+2)=LX(3,j+1); 
        LX2(3,k+2)=LX(3,j+2); 
        LX2(4,k+2)=LX(4,j+2); 
   
        LY2(1,k+2)=LY(1,j+1); 
        LY2(2,k+2)=LY(3,j+1); 
        LY2(3,k+2)=LY(3,j+2); 
        LY2(4,k+2)=LY(4,j+2);

        iflocked(k+2)=1;

%       Delete K and K+1th element
        LX2(:,[k k+1]) = [];
        LY2(:,[k k+1]) = [];
        
        iflocked([k k+1]) = [];

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

%    l1=length(LX2);
%    cmap=jet(l1);
%    for j=1:l1
%      xt=LX2(:,j);
%      yt=LY2(:,j);
%      figure(fig2);
%      fill(xt,yt,cmap(j,:)); hold on
%    end
  end  

end   % end function

%---------------------------------------------------------------------- 





