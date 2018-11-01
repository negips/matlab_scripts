function ind = FaceIndex3D(cf)

  if cf==1
    ind=[1 2 6 5];
  elseif cf==2
    ind=[2 3 7 6];
  elseif cf==3
    ind=[3 4 8 7];
  elseif cf==4
    ind=[4 1 5 8];
  elseif cf==5
    ind=[1 2 3 4];
  elseif cf==6
    ind=[5 6 7 8];
  else
    disp(['Unknown face index value: ', num2str(cf)])
    ind = [];
  end  

end   % function
%---------------------------------------------------------------------- 
