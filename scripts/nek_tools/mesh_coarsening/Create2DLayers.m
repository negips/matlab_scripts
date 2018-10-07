function [NewE, NewX, NewY, NewBC, NewCEl, NewCoF, NewET, iflocked]=Create2DLayers(NewE,NewX,NewY,NewBC,NewCEl,NewCoF,NewET,cr_layer,ifc,if3skip,fig2,ifplot)

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
    LET=NewET{i};
 
    LE2=LE; 
    LX2=LX;
    LY2=LY;
    LBC2=LBC;
    LCEl2=LCEl;
    LCoF2=LCoF;
    LET2=LET;
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

        LET2{k-1}    = 'e3';

        iflocked(k-1)=1;
        if k>2
          iflocked(k-2)=1;
        end  

%       Enlarge K+2th element
        if j<=l1-2
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

          LET2{k+2}    = 'e1';

          iflocked(k+2)=1;
        end  

%       Delete K and K+1th element
        if j<=l1-2
          LE2([k k+1])   = [];
          LX2(:,[k k+1]) = [];
          LY2(:,[k k+1]) = [];
          LBC2([k k+1])  = [];
          LCEl2(:,[k k+1])= [];
          LCoF2(:,[k k+1])= [];
          LET2([k k+1])   = [];
    
          iflocked([k k+1]) = [];
        else
          LE2(k)   = [];
          LX2(:,k) = [];
          LY2(:,k) = [];
          LBC2(k)  = [];
          LCEl2(:,k)= [];
          LCoF2(:,k)= [];
          LET2(k)   = [];

          iflocked(k) = [];
        end

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
        if j<=l1-2
          LCEl2(3,k-1) = LCEl(3,j+1);
        else
          LCEl2(3,k-1) = LCEl(3,j);
        end  

        LCEl2(4,k-1) = LCEl(4,j-1);

        LCoF2(1,k-1) = LCoF(1,j-1);
        LCoF2(2,k-1) = LCoF(2,j-1);
        if j<=l1-2
          LCoF2(3,k-1) = LCoF(3,j+1);
        else  
          LCoF2(3,k-1) = LCoF(3,j);
        end  

        LCoF2(4,k-1) = LCoF(4,j-1);

        LET2{k-1}    = 's';   % no change


%       Enlarge K+2th element
        if j<=l1-2
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

          LET2{k+2}    = 's';   % no change
        end

%       Delete K and K+1th element
        if j<=l1-2
          LE2([k k+1])   = [];
          LX2(:,[k k+1]) = [];
          LY2(:,[k k+1]) = [];
          LBC2([k k+1])  = [];
          LCEl2(:,[k k+1])= [];
          LCoF2(:,[k k+1])= [];
          LET2([k k+1])   = [];
    
        else
          LE2(k)   = [];
          LX2(:,k) = [];
          LY2(:,k) = [];
          LBC2(k)  = [];
          LCEl2(:,k)= [];
          LCoF2(:,k)= [];
          LET2(k)   = [];

        end

      end

      if j<=l1-2
        k=k-2;
      else
        k=k-1;
      end  

    end 

    NewE{i}=LE2;
    NewX{i}=LX2;
    NewY{i}=LY2;
    NewBC{i}=LBC2;
    NewCEl{i}=LCEl2;
    NewCoF{i}=LCoF2;
    NewET{i}=LET2;

    if (ifplot)
      l1=length(LX2);
      cmap=jet(l1);
      for j=1:l1
        xt=LX2(:,j);
        yt=LY2(:,j);
        figure(fig2);
        fill(xt,yt,cmap(j,:)); hold on
      end
    end  
  end  

end   % end function

%---------------------------------------------------------------------- 


