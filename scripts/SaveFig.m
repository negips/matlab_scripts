function SaveFig(h, filename, destn, col)
% Save and move the figure


if ~isempty(strfind(lower(filename), '.eps'))
  if col
       print( h, '-r150', '-depsc2', filename)
  else
       print( h, '-r150', '-deps2', filename)
  end
  movefile([filename],destn)

elseif ~isempty(strfind(lower(filename), '.pdf'))
  print( h, '-r150', '-dpdf', filename)
  movefile([filename],destn)

elseif ~isempty(strfind(lower(filename), '.png'))
  print( h, '-r150', '-dpng', filename)
  movefile([filename],destn)

else
  if col
       print( h, '-r150', '-depsc2', filename)
  else
       print( h, '-r150', '-deps2', filename)
  end

  movefile([filename '.eps'],destn)
end


end

