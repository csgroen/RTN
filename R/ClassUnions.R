##definitions of class unions
setClassUnion("num_Or_int_Or_null",c("numeric","integer","NULL"))
setClassUnion("num_Or_int",c("numeric","integer"))
setClassUnion("char_Or_null", c("character", "NULL"))
setClassUnion("num_Or_null", c("numeric", "NULL"))
