
% student-t distributed noise
function samples = trandn(m,n,dof)
    variance = dof/(dof-2);
    coef = 1/sqrt(variance);
    
    samples = coef * trnd(dof,m,n);     % "samples" has zero mean and unit variance
end