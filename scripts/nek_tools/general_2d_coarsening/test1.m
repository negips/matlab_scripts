%     Trying out a general script to do 2D coarsening
%    

addpath '/scratch/negi/git_repos/matlabscripts/scripts/nek_tools/mesh_coarsening/'

casename = 'vlarge1';

%rea = Nek_ReadRea(casename)

steln = 2000;        % Element to start coarsening
stedg = 1;           % which edge to coarsen

cbc = rea.mesh.cbc;
% Find neighbouring elements that need modification

el3no = cbc(stedg,steln).connectsto;        % Top element
if el3no>0
  el3 = cbc(:,el3no);
else
  el3 = [];
end  

edgl = CorrectEdgeNo(stedg-1);              % Left Edge no
el1no = cbc(edgl,steln).connectsto;         % Left element
if el1no>0
  el1 = cbc(:,el1no);
else
  el1 = [];
end  

edgr = CorrectEdgeNo(stedg+1);              % Left Edge no
el7no = cbc(edgr,steln).connectsto;         % Right element
if el7no>0
  el7 = cbc(:,el7no);
else
  el7 = [];
end  







function edge = CorrectEdgeNo(edge)

   if edge>4
     edge = edge-4;
   elseif edge<1
     edge = edge+4;
   end

end 

