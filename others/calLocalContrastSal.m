function cocomSal = calLocalContrastSal(W_lab,bgcomSal,disSup,sal_lab,spnum,adjc,theta,num_vals,x_vals,y_vals)

disSup = normalize( disSup );
sal_lab = sal_lab + eye(spnum);
DB = W_lab.*adjc;
[val,index] = max(DB,[],2);
ind = find(val > mean(val));
index = index(ind);
conLab = sal_lab(ind,:);
simLab = conLab.*exp( -theta*disSup(ind,:) );

salSup = simLab.*(ones(length(ind),1)*num_vals');
Isum = sum(salSup,2);
x_valMat = ones(length(ind),1)*x_vals'; y_valMat = ones(length(ind),1)*y_vals';
clear DB conLab

cocomSal = zeros(spnum,1);
% É¢²¼³Ì¶È
sal = bgcomSal; sal(bgcomSal < 0.5) = 0;
salCom = sal.*num_vals;
IsumCom = sum(salCom);
x0 = sum(salCom.*x_vals)/IsumCom;
y0 = sum(salCom.*y_vals)/IsumCom;
coherence = sum(salSup.*sqrt((x_valMat-x0).^2 + (y_valMat - y0).^2),2)./Isum;

cocomSal(ind) = val(ind).*(1 - normalize(coherence));
cocomSal(ind) = cocomSal(ind).*( cocomSal(ind) > cocomSal(index) );

cocomSal = max((cocomSal(ind)*ones(1,spnum)).*simLab,[],1)';
cocomSal(cocomSal < mean(cocomSal)) = mean(cocomSal);
cocomSal = normalize( cocomSal );
clear sal salCom IsumCom x0 y0 coherence simLab salSup x_valMat y_valMat val index Isum