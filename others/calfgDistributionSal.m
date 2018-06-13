function comSal = calfgDistributionSal(salSup,Isum,x_vals,y_vals,m,n,spnum)

dist = calDistribution(salSup,Isum,x_vals,y_vals,m,n,spnum,spnum);

comSal = dist';
comSal(comSal > mean(comSal)) = mean(comSal);
comSal = 1 - normalize(comSal);
comSal = comSal';
clear dist coherence centric
