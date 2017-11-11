function SaveFig(h, filename, destn, col)
% Save and move the figure

reso='-r300';


if ~isempty(strfind(lower(filename), '.eps'))
  if col
       print( h, reso, '-depsc2', filename)
  else
       print( h, reso, '-deps2', filename)
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
       print( h, reso, '-depsc2', filename)
  else
       print( h, reso, '-deps2', filename)
  end

  movefile([filename '.eps'],destn)
end


end

