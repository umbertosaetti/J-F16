function fix(a)
    if a>0
        return floor(a)
    elseif a<0
        return ceil(a)
    else
        return a
    end
end 
