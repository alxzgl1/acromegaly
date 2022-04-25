function res = iRoman(number)
number = lower(number);
res = 0;
for i = 1:length(number)
    if number(i) == 'i' && i < length(number) && contains(number(i+1:length(number)), {'v','x'})
        res = res - 1;
    elseif number(i) == 'i'
        res = res + 1;
    elseif number(i) == 'v' && i < length(number) && contains(number(i+1:length(number)), {'x'})
        res = res - 5;
    elseif number(i) == 'v'
        res = res + 5;
    elseif number(i) == 'x'
        res = res + 10;
    end
end
end