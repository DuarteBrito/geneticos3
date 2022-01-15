function [a,b] = sig(infe, sup)
    inte = sup-infe;
    infe = infe-2*inte;
    b = -8/(sup-infe);
    a = -((sup-infe)/2+infe)*b;
end

