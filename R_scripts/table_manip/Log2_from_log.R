########## getting log '2' FC from logFC #############



test_S8vsD8 <- All_S8vsD8

test_S8vsI8 <- All_S8vsI8

test_I8vsD8 <- All_I8vsD8

test_I8vsD8$log2FC <- (test_I8vsD8$logFC/log(2))
