function SaveFig(h, filename, destn, col)
% Save and move the figure

reso='-r100';


if ~isempty(strfind(lower(filename), '.eps'))
  if col
       print( h, reso, '-depsc2', filename)
  else
       print( h, reso, '-deps2', filename)
  end
  if ~isempty(destn) 
    movefile([filename],destn)
  end  

elseif ~isempty(strfind(lower(filename), '.pdf'))
  print( h, reso, '-dpdf', filename)
  if ~isempty(destn) 
    movefile([filename],destn)
  end  

elseif ~isempty(strfind(lower(filename), '.png'))
  print( h, reso, '-dpng', filename)
  if ~isempty(destn) 
    movefile([filename],destn)
  end  

else
  if col
       print( h, reso, '-depsc2', filename)
  else
       print( h, reso, '-deps2', filename)
  end

  if ~isempty(destn) 
    movefile([filename '.eps'],destn)
  end  
 
end


end

