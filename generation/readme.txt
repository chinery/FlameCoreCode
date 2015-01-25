To create a new flame animation, run makeanimation.m

The script loads in some training data, stored in ../data/generation

There is some horrible hackery going on in this code to make my life easier coding it. The same routine is performed for several parameters, so instead of doing this:
param1 = ...
param2 = ...
param3 = ...
...

it does
params = {'param1', 'param2', 'param3',...}
for i = 1:length(params)
	eval(sprintf('%s = ...',params{i});
end

which is not the correct way of handling these sorts of structures, but it works.

Explanation of parameters:
parameters.animation.totalitr  				-- how many animations to create with the same settings
parameters.animation.anilength 				-- animation length in frames
parameters.animation.start 					-- starting frame number
parameters.animation.dirstring 				-- where to store the results, include a %i for the loop over all animations
parameters.nearest.num 						-- number of nearest points, from training data, to consider when calculating where to go next
parameters.nearest.distance 				-- but don't go further than this number of standard deviations away
parameters.nearest.exactdistcorrection		-- if there are any points of distance 0, give them some distance so probability isn't Inf%
parameters.gaussian.num 					-- generate this many points around the target during randomisation
parameters.gaussian.scale					-- scale the gaussian during randomisation, higher number = further away from training
											   where 1 is the average distance between points, so probably should be 0-0.5
modifier.param 								-- some parameters handle the randomisation better than others, so can be manually adjusted