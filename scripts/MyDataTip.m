function outText = MyDataTip(obj, event)
  pos = get(event, 'Position');
  if length(pos) == 3
    outText = {sprintf('X: %.5f', pos(1), ...
               sprintf('Y: %.5f', pos(2), ...
               sprintf('Z: %.5f', pos(3)};
  else
    outText = {sprintf('X: %.5f', pos(1), ...
               sprintf('Y: %.5f', pos(2)};
  end

return
