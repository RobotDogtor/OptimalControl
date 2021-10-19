function checkIfV123Failed(V1,V2,V3,epsilon)
V1fail = 0;
V2fail = 0;
V3fail = 0;
for i = 1:(length(V1)-1)
    if (V1(i+1)-V1(i))>epsilon
        V1fail = V1fail+ 1;
    end
end
for i = 1:(length(V2)-1)
    if (V2(i+1)-V2(i))>epsilon
        V2fail = V2fail+1;
    end
end
for i = 1:(length(V3)-1)
    if (V3(i+1)-V3(i))>epsilon
        V3fail = V3fail+1;
    end
end
if V1fail>1
     disp(['Lyapunov Function 1 failed: ' num2str(V1fail) '/' num2str(length(V1))]);
end
if V2fail>1
     disp(['Lyapunov Function 2 failed: ' num2str(V2fail) '/' num2str(length(V2))]);
end
if V3fail>1
     disp(['Lyapunov Function 3 failed: ' num2str(V3fail) '/' num2str(length(V3))]);
end
if (V1fail+V2fail+V3fail)==0
    disp(['All 3 Lyapunov Derivatives are negative, to epsilon = ' num2str(epsilon)])
end

eV = abs(V1-V2);
disp(['avarage error between V1 and V2: ' num2str(mean(eV))]);


end