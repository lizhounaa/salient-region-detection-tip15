function comSal = calfgCompactnessSal(salSup,Isum,disSup,x_vals,y_vals,weight,m,n,spnum)

% xo = n/2; yo = m/2;

% disCen = sqrt((ones(spnum,1)*x_vals'-xo).^2 + (ones(spnum,1)*y_vals'-yo).^2);
coh = sum(salSup.*disSup,2)./Isum;
% cen = sum(salSup.*disCen,2)./Isum;

bgDis = min([y_vals-1 n-x_vals m-y_vals x_vals-1],[],2);
bgDis = max(bgDis) - bgDis;
cen = sum(salSup.*(ones(spnum,1)*bgDis'),2)./Isum;

clear xo yo salSup Isum disCen

comSal = coh'+cen';%normalize(coh'+cen');
comSal(comSal > mean(comSal)) = mean(comSal);
comSal = 1 - normalize(comSal);
comSal = normalize( sum(weight.*padarray(comSal,[spnum-1 0],'replicate','post'),2) );
clear coh cen
