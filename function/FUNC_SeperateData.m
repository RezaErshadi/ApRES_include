function [Signal] = FUNC_SeperateData(Data)

for i = 1:size(Data,2)
    or(1,i) = string(Data{2,i})+string(Data{3,i})+string(Data{4,i});
end

unqor = unique(or,'stable');

fprintf("----------\nSignal Report:\n")

for i = 1:length(unqor)
    temp2 = find(or == unqor(i));
    fprintf("%s: %d data\n",unqor(i),length(temp2))
    for j = 1:length(temp2)
        ix = temp2(j);
        eval("Signal."+unqor(i)+"(:,"+string(j)+") = Data{end,"+string(ix)+"};"); 
    end
end
