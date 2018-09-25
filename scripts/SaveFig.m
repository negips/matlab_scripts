function SaveFig(h, filename, destn, col)
% Save and move the figure

reso='-r100';

renderer='-opengl';
%renderer='-painters';


if ~isempty(strfind(lower(filename), '.eps'))
  if col
       print( h, renderer, reso, '-depsc2', '-tiff', filename)
  else
       print( h, renderer, reso, '-deps2', '-tiff', filename)
  end
  movefile([filename],destn)

elseif ~isempty(strfind(lower(filename), '.pdf'))
  print( h, reso, '-dpdf', filename)
  movefile([filename],destn)

elseif ~isempty(strfind(lower(filename), '.png'))
  print( h, reso, '-dpng', filename)
  movefile([filename],destn)

else
  if col
       print( h, renderer, reso, '-depsc2', filename)
  else
       print( h, renderer, reso, '-deps2', filename)
  end

  movefile([filename '.eps'],destn)
end


end

