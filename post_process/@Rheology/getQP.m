% QP = GETQP(obj) get time averaged deviatric pressure q/p
%
% For each time, read the cell stresses and compute the principle stresses eigs.
% Then find the deviatric pressure q/p, where q = (eigs(1)-eigs(2))/2
 % and p = (eigs(1)+eigs(1))/2

function qp = getQP(obj)

    eigs = eig(obj.stress);
    qp = (eigs(1)-eigs(2))/(eigs(2)+eigs(1));

end
