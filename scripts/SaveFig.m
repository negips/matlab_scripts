function SaveFig(h, filename, destn, col)
% Save and move the figure

reso='-r150';

%renderer='-opengl';
renderer='-painters';


if ~isempty(strfind(lower(filename), '.eps'))
  if col
       print( h, renderer, reso, '-depsc2', '-tiff', filename)
  else
       print( h, renderer, reso, '-deps2', '-tiff', filename)
  end
  if ~isempty(destn) 
    movefile([filename],destn)
  end  

elseif ~isempty(strfind(lower(filename), '.pdf'))
  print(h, reso, '-dpdf', '-bestfit', filename)
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
       print( h, renderer, reso, '-depsc2', filename)
  else
       print( h, renderer, reso, '-deps2', filename)
  end

  if ~isempty(destn) 
    movefile([filename '.eps'],destn)
  end  
 
end

disp(['Figure Saved: ', destn, ' ', filename])

end

