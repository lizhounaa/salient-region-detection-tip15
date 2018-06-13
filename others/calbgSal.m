function bgSal = calbgSal(comSal,sal_lab,x_vals,y_vals,weight,m,n,spnum)

% ±³¾°Ë÷Òý
bgIndex = find(comSal < mean(comSal)); 
sal = 1 - comSal(bgIndex);
bgdisSal_lab = sal_lab(bgIndex,:);
bgX = x_vals(bgIndex); bgY = y_vals(bgIndex);
[~,index_dis] = min([bgY-1 n-bgX m-bgY bgX-1],[],2);
bgNum_dis = unique(index_dis);
bgSal = ones(1,spnum);
for i = 1:length(bgNum_dis)
    cs = sal(index_dis==bgNum_dis(i));
    [val, ind] = min(1-bgdisSal_lab(index_dis==bgNum_dis(i),:),[],1);
    if length(cs) == 1
        tempSal = normalize( val );
    else
        tempSal = normalize( val.*cs(ind)' );
    end
    bgSal = bgSal.*tempSal;
    clear tempSal cs val ind
end
bgSal = normalize( sum(weight.*padarray(bgSal,[spnum-1 0],'replicate','post'),2) );
clear bgIndex bgX bgY sal

