% First attempts at coarsening the mesh

clear
clc
close all

load saab_wing2d.mat

skiplayers = 5;         % No of layers to skip when coarsening
lmax = 0.1;             % Maximum length of a side. If the length is larger. Don't coarsen
ARcut = 2.5;            % Coarsen if Aspect ratio is larger than this.


% invert numbering and make new arrays/cells

for i=1:nlayers
  j=nlayers-i+1;
  NewEl{i}=LayersEl{j};
  NewFV{i}=LayersFopV{j};
  NewFO{i}=LayersFopO{j};
end


% Test coarsening algorithm 1
% Coarsen Layer by Layer
% Define aspect ratio as 'O' face lengths to 'V' face lengths
layer_start = skiplayers+1;
for i=layer_start:layer_start  %nlayers

  Lel=NewEl{i};
  LV =NewFV{i};
  LO =NewFO{i};  
  l1 =length(Lel);

  cmap = jet(l1); 
  for j=1:l1

    e = Lel(j);
        
    % Index of 'O' elements
    k   = LO(j);
    k1 = k;
    k2 = k+1;
    if k2>4
      k2=k2-4;
    end
    ko1 = [k1; k2];
    dx  = diff(rea.mesh.xc(ko1,e));
    dy  = diff(rea.mesh.yc(ko1,e));
    dl1 = sqrt(dx^2 + dy^2);

    k   = LO(j);
    k   = k+2;
    if (k>4)
      k=k-4;
    end  
    k1 = k;
    k2 = k+1;
    if k2>4
      k2=k2-4;
    end
    ko2 = [k1; k2];
    dx  = diff(rea.mesh.xc(ko2,e));
    dy  = diff(rea.mesh.yc(ko2,e));
    dl2 = sqrt(dx^2 + dy^2);
    dlo = mean([dl1 dl2]);  

    % Index of 'v' elements
    k   = LV(j);
    k1 = k;
    k2 = k+1;
    if k2>4
      k2=k2-4;
    end
    kv1 = [k1; k2];
    dx  = diff(rea.mesh.xc(kv1,e));
    dy  = diff(rea.mesh.yc(kv1,e));
    dl1 = sqrt(dx^2 + dy^2);

    k   = LV(j);
    k   = k+2;
    if (k>4)
      k=k-4;
    end  
    k1 = k;
    k2 = k+1;
    if k2>4
      k2=k2-4;
    end
    kv2 = [k1; k2];
    dx  = diff(rea.mesh.xc(kv2,e));
    dy  = diff(rea.mesh.yc(kv2,e));
    dl2 = sqrt(dx^2 + dy^2);
    dlv = mean([dl1 dl2]);

    l_ar(j) = dlo/dlv;

    xt = rea.mesh.xc(:,e);
    yt = rea.mesh.yc(:,e);
    fill(xt,yt,cmap(j,:)); hold on

  end 

  



end



