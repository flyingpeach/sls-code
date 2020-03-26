function bool = is_instance(obj, classStr)
    bool = strcmp(class(obj), classStr);
end