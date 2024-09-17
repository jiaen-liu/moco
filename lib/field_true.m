function y=field_true(struct,field)
    y=(isfield(struct,field) && getfield(struct,field));
end