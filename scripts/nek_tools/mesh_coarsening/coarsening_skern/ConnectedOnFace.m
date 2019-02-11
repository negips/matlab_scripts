function NewCoF = ConnectedOnFace(NewE,NewCEl,NewMeshC)

  nl=length(NewMeshC);              % no of layers
  
  selfwrap  = 1;            % This layer wrap around on itself  
  for il=1:nl
        
    newcof = []; 
    [r nel] = size(NewCEl{il});
    nfaces=r;
    foundwall = 0;            % Needed for first C type layer. onface value changes before/after the airfoil

    for i=1:nel
      f=int32(zeros(nfaces,1));
      for j=1:nfaces

        if j==2 && NewCEl{il}(j,i)==0
          foundwall=1;
        end 

        if NewCEl{il}(j,i)==0   % Boundary element
          f(j)=0;
        else
          f2=j+2;
          if f2>4
            f2=f2-4;
          end
          f(j)=f2;

          if (j==2 && il==selfwrap)
            f(j)=2;
            continue
          end

          if j==4 || j==2 
%           find layer type of connecting element
            cel = NewCEl{il}(j,i);
            for il2=1:nl
              ind=find(NewE{il2}==cel);
              if ~isempty(ind)
                onlayer=il2;
                break
              end
            end
            onLC = NewMeshC(onlayer);

            if j==2 && NewMeshC(il)==0 && onLC==1  % This layer is not C type, Connecting layer is C type
               f(j)=2;
            elseif j==2 && NewMeshC(il)==1 && onLC==0  % This layer is C type, Connecting layer is not C type
               if foundwall
                 f(j)=2;
               else
                 f(j)=4;
               end
            end   
          end     % j==4

        end     % NewCEl{il}(j,i)==0
    
      end       % j=1:nfaces
      newcof = [newcof f]; 
    end         % i=1:nel
    NewCoF{il}=newcof;
  end           % il   

end   % function  
%---------------------------------------------------------------------- 


