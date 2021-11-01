function value = phi(normalized_innovation, a,b,c)
    m = length(normalized_innovation);
    value = normalized_innovation;
    
    for j = 1:m
		if normalized_innovation(j) >= 0
			if normalized_innovation(j) <= a
				%warning('here 1')
				value(j) = c*tan(c*normalized_innovation(j)/2);
			elseif normalized_innovation(j) >= b
				%warning('here 2')
				value(j) = b;
			end
		else
			if normalized_innovation(j) >= -a
				%warning('here 3')
				value(j) = -c*tan(-c*normalized_innovation(j)/2);
			elseif normalized_innovation(j) <= -b
				%warning('here 4')
				value(j) = -b;
			end
		end
    end
end