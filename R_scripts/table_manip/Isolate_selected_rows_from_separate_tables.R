#### This line allows you to create a new table from selecting entire rows in one table based on values of a column in a separate table ######

###### the "values_1" represent genes IDS (for example) which are the same in both tables however "seperate_table" has only the ones we want
######

new_table <- old_table[is.element(old_table$values_1, separate_table$values_1),]