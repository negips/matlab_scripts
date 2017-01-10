function SaveFig(h, filename, destn, col)
% Save and move the figure

if col

     print( h, '-r600', '-depsc', filename)
else
     print( h, '-r600', '-deps', filename)
end

movefile([filename '.eps'],destn)

end

