function [a,b] = remove_multiple_entries(a,b)

num_a = length(a);
i=1;
j=1;

while i <= num_a
    j=i;
    while j <= num_a
        if abs(a(i)-a(j)) < 1e-6  && abs(b(i)-b(j))<1e-6 && i ~=j
            a(j) = [];
            b(j) = [];
            num_a = num_a - 1;
        end
        j=j+1;
    end
    i=i+1;
end