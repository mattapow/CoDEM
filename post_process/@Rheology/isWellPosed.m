function isWP = isWellPosed(f, I)
	% isWP = isWellPosed(f, I) 
	% Wether the given mu(I) fit (f) is linearly well posed at I
	% 
	% See equation 2.16 for the incompressible governing equations:
	% Well-posed continuum equations for granular  ow with compressibility and μ(I)-rheology
	% http://dx.doi.org/10.1098/rspa.2016.0846

	assert(strcmp(formula(f), 'a+(b-a)/(c/x+1)'))

	mu = f(I);
	dMu = differentiate(f, I);
	v = I/mu * dMu; 

	V = 4*v^2 - 4*v + mu^2 * (1 - v/2)^2;

	if V <= 0
		isWP = 1;
	else
		isWP = 0;
	end	 

end