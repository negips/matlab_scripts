% First attempts at coarsening the mesh

clear
clc
close all


%           Element arrangement
%
%
%                 f2
%           x3-----------x2      
%           |            |
%           |            |
%         f3|            |f1 'O  '
%           |            |
%           |            |
%           x4-----------x1
%                 f4
%               'v  '


load saab_wing2d.mat
%load saab750k.mat

skiplayers = 2;         % Need to skip first layer since its smaller than the others

disp(['Total Number of Elements: ', num2str(rea.mesh.nelg)])

for i=1:nlayers
  j=nlayers-i+1;
  OldE{i}=LayerE{j};
  OldX{i}=LayerX{j};
  OldY{i}=LayerY{j};
  OldBC{i}=LayerBC{j};
  OldCEl{i}=LayerCEl{j};

  NewE{i}=OldE{i};
  NewX{i}=OldX{i};
  NewY{i}=OldY{i};
  NewBC{i}=OldBC{i};
  NewCEl{i}=OldCEl{i};
  NewCoF{i}=ConnectedOnFace(NewCEl{i});

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
  LCEl=NewCEl{i}; 

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

% Coarsen layer in consecutive pairs
  figure(2);    
  [LX,LY,ifc,NewX,NewY,NewCEl,NewCoF] = CoarsenLayer(LE,LX,LY,ifc,NewX,NewY,NewCEl,NewCoF,i,fig2);
% Modify all subsequent Layers
  if3skip=1;  
  [NewE, NewX, NewY, NewBC, NewCEl, NewCoF, iflocked]=CreateNewLayers(NewE,NewX,NewY,NewBC,NewCEl,NewCoF,i,ifc,if3skip,fig2);
  nels_layer = length(NewX{i});
  new_nelg = new_nelg+nels_layer;
end
disp(['Total Number of Elements: ', num2str(new_nelg)])

% Not done for multiple definitions right now
curvedef= 'mv ';
Mesh2 = ReOrderElements(NewE,NewX,NewY,NewBC,NewCEl,NewCoF,rea.mesh,curvedef); 

%nz0=5;
%Lz=1.0;
%[NewE3, NewX1, NewY1, NewX2, NewY2, NewBC3, NewCEl3] = Generate3D(NewE,NewX,NewY,NewBC,NewCEl,nlayers,nz0,Lz);
%---------------------------------------------------------------------- 

function newcof = ConnectedOnFace(NewCEl)

  [r nel] = size(NewCEl);
 
  newcof = []; 
  nfaces=r;
  for i=1:nel
    f=int32(zeros(nfaces,1));
    for j=1:nfaces
      if NewCEl(j,i)==0   % Boundary element
        f(j)=0;
      else
        f2=j+2;
        if f2>4
          f2=f2-4;
        end
        f(j)=f2;
      end     % if
  
    end       % j=1:nfaces
    newcof = [newcof f]; 
  
  end         % i=1:nel

end   % function  
%---------------------------------------------------------------------- 

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

function [LX,LY,ifc,NewX,NewY,NewCEl,NewCoF] = CoarsenLayer(LE,LX,LY,ifc,NewX,NewY,NewCEl,NewCoF,i,fig2)

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
      figure(fig2)
      xt = LX(:,j);
      yt = LY(:,j);
      figure(2)
      fill(xt,yt,cmap(j,:)); hold on
      
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

%   Coarsen Element j+1
    LX(1,j+1)=LX(4,j);
    LY(1,j+1)=LY(4,j);
%   Connecting Element no changes  
    NewCEl{i}(4,j+1)=NewCEl{i}(4,j+2); 
%   onFace of connecting element changes
    NewCoF{i}(4,j+1)=1;  

    figure(fig2)
    xt = LX(:,j);
    yt = LY(:,j);
    figure(2)
    fill(xt,yt,cmap(j,:)); hold on

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

function [NewE, NewX, NewY, NewBC, NewCEl, NewCoF, iflocked]=CreateNewLayers(NewE,NewX,NewY,NewBC,NewCEl,NewCoF,cr_layer,ifc,if3skip,fig2)

  iflocked = [];

  nlayers=length(NewX);
  il=0;
  for i=cr_layer+1:nlayers
    il=il+1;  
    LE=NewE{i};
    LX=NewX{i};
    LY=NewY{i};
    LBC=NewBC{i};
    LCEl=NewCEl{i};
    LCoF=NewCoF{i};
 
    LE2=LE; 
    LX2=LX;
    LY2=LY;
    LBC2=LBC;
    LCEl2=LCEl;
    LCoF2=LCoF;
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

        LBC2{k-1}(1,:) = LBC{j-1}(1,:);
        LBC2{k-1}(2,:) = LBC{j-1}(2,:);
        LBC2{k-1}(3,:) = LBC{j}(2,:);
        LBC2{k-1}(4,:) = LBC{j-1}(4,:);

        LCEl2(1,k-1) = LCEl(1,j-1);
        LCEl2(2,k-1) = LCEl(2,j-1);
        LCEl2(3,k-1) = LCEl(2,j);
        LCEl2(4,k-1) = LCEl(4,j-1);

        LCoF2(1,k-1) = LCoF(1,j-1);
        LCoF2(2,k-1) = LCoF(2,j-1);
        LCoF2(3,k-1) = 4;
        LCoF2(4,k-1) = LCoF(4,j-1);

        iflocked(k-1)=1;
        if k>2
          iflocked(k-2)=1;
        end  

%       Enlarge K+2th element        
        LX2(1,k+2)=LX(1,j+1); 
        LX2(2,k+2)=LX(3,j+1); 
        LX2(3,k+2)=LX(3,j+2); 
        LX2(4,k+2)=LX(4,j+2); 
   
        LY2(1,k+2)=LY(1,j+1); 
        LY2(2,k+2)=LY(3,j+1); 
        LY2(3,k+2)=LY(3,j+2); 
        LY2(4,k+2)=LY(4,j+2);

        LBC2{k+2}(1,:) = LBC{j+1}(2,:);
        LBC2{k+2}(2,:) = LBC{j+2}(2,:);
        LBC2{k+2}(3,:) = LBC{j+2}(3,:);
        LBC2{k+2}(4,:) = LBC{j+2}(4,:);

        LCEl2(1,k+2) = LCEl(2,j+1);
        LCEl2(2,k+2) = LCEl(2,j+2);
        LCEl2(3,k+2) = LCEl(3,j+2);
        LCEl2(4,k+2) = LCEl(4,j+2);

        LCoF2(1,k+2) = 4;
        LCoF2(2,k+2) = LCoF(2,j+2);
        LCoF2(3,k+2) = LCoF(3,j+2);
        LCoF2(4,k+2) = LCoF(4,j+2);


        iflocked(k+2)=1;

%       Delete K and K+1th element
        LE2([k k+1])   = [];
        LX2(:,[k k+1]) = [];
        LY2(:,[k k+1]) = [];
        LBC2([k k+1])  = [];
        LCEl2(:,[k k+1])= [];
        LCoF2(:,[k k+1])= [];
     
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

        LBC2{k-1}(1,:) = LBC{j-1}(1,:);
        LBC2{k-1}(2,:) = LBC{j-1}(2,:);
        LBC2{k-1}(3,:) = LBC{j}(3,:);
        LBC2{k-1}(4,:) = LBC{j-1}(4,:);

        LCEl2(1,k-1) = LCEl(1,j-1);
        LCEl2(2,k-1) = LCEl(2,j-1);
        LCEl2(3,k-1) = LCEl(3,j+1);
        LCEl2(4,k-1) = LCEl(4,j-1);

        LCoF2(1,k-1) = LCoF(1,j-1);
        LCoF2(2,k-1) = LCoF(2,j-1);
        LCoF2(3,k-1) = LCoF(3,j);
        LCoF2(4,k-1) = LCoF(4,j-1);

%       Enlarg K+2th element        
        LX2(1,k+2)=LX(1,j+1); 
        LX2(2,k+2)=LX(2,j+1); 
        LX2(3,k+2)=LX(3,j+2); 
        LX2(4,k+2)=LX(4,j+2); 
   
        LY2(1,k+2)=LY(1,j+1); 
        LY2(2,k+2)=LY(2,j+1); 
        LY2(3,k+2)=LY(3,j+2); 
        LY2(4,k+2)=LY(4,j+2);

        LBC2{k+2}(1,:) = LBC{j+1}(1,:);
        LBC2{k+2}(2,:) = LBC{j+2}(2,:);
        LBC2{k+2}(3,:) = LBC{j+2}(3,:);
        LBC2{k+2}(4,:) = LBC{j+2}(4,:);

        LCEl2(1,k+2) = LCEl(1,j);
        LCEl2(2,k+2) = LCEl(2,j+2);
        LCEl2(3,k+2) = LCEl(3,j+2);
        LCEl2(4,k+2) = LCEl(4,j+2);

        LCoF2(1,k+2) = LCoF(1,j+1);
        LCoF2(2,k+2) = LCoF(2,j+2);
        LCoF2(3,k+2) = LCoF(3,j+2);
        LCoF2(4,k+2) = LCoF(4,j+2);

%       Delete K and K+1th element
        LE2([k k+1])   = [];
        LX2(:,[k k+1]) = [];
        LY2(:,[k k+1]) = [];
        LBC2([k k+1])  = [];
        LCEl2(:,[k k+1])= [];
        LCoF2(:,[k k+1])= [];

      end

      k=k-2;

    end 

    NewE{i}=LE2;
    NewX{i}=LX2;
    NewY{i}=LY2;
    NewBC{i}=LBC2;
    NewCEl{i}=LCEl2;
    NewCoF{i}=LCoF2;

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





