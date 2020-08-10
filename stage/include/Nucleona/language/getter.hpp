#define GETTER(return_type, getter_sig, statement ) \
    return_type getter_sig \
        statement \
    const return_type getter_sig const \
        statement

#define CGETTER(return_type, getter_sig, statement ) \
    const return_type getter_sig const \
        statement

#define T_GETTER(tpl, return_type, getter_sig, statement ) \
    tpl \
    return_type getter_sig \
        statement \
    tpl \
    const return_type getter_sig const \
        statement
#define T_CGETTER(tpl, return_type, getter_sig, statement ) \
    tpl \
    const return_type getter_sig const \
        statement

