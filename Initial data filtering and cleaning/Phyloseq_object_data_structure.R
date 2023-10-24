#code to change metadata info of pseq object


#change data type ----

# Convert the column 'YourColumn' to a character data type

pseq.core@sam_data$YourColumn <- as.character(pseq.core@sam_data$YourColumn)

pseq.core@sam_data

