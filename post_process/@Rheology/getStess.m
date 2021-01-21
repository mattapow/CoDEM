% OBJ = GETSTRESS(OBJ) Get time averaged cell stress

function obj = set.stess(obj)

if ~any(isnan(obj.stress))
  return
end

for filenum = obj.tEnd:-1:obj.nPass
    % Read in cell stress
    fp = strcat(obj.dirPath, '/cell/cell_', string(filenum));
    cell = Cell();
    formatspec = '%*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %*f %*f';
    cell = readCell(cell, fp, formatspec);
    stress(filenum-obj.nPass+1, :) = cell.stress;
end

assert(~any(stress==0), 'Zeros stresses detected')
assert(~any(isnan(stress)), 'Zeros stresses detected')

obj.stress = mean(stress, 1);

end
