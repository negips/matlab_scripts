% Keeping the script separate for single layer coarsening


l2=length(ifc);
eltmp = [];
fvtmp = [];
fotmp = [];
xct   = [];
yct   = [];
for kk=1:l2
  if (~ifc(kk))
    eltmp = [eltmp Lel(kk)];
    fvtmp = [fvtmp LV(kk)];
    fotmp = [fotmp LO(kk)];
    xct   = [rea.mesh.xc(:,Lel(kk))];
    yct   = [rea.mesh.yc(:,Lel(kk))];
    continue
  end

  % face opposite 'O'
  f = LO(kk);
  % Not sure how orientations are specified.
  % So I just test the vertices of this face with the vertices of the next layer elements.
  f1 = f;
  f2 = f+1;
  if (f2>4)
    f2=f2-4;
  end
  

end  

  

    
