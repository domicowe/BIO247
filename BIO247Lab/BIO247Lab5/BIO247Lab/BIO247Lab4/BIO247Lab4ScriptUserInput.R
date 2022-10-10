#Script for Removing Junk LEEs: USER INPUT
#Whitney Domico BIO247 Lab 4

#user input file path
file <- readline(prompt="Enter file path: ")

#upload file as dataframe
filename <- read.csv(file,header=TRUE)

#create vector of paper IDs
LEEs <- filename$Paper.ID

#remove "junk" characters
filename$Paper.ID <- gsub("\\]", "", filename$Paper.ID)
filename$Paper.ID <- gsub("\\[", "", filename$Paper.ID)
filename$Paper.ID <- gsub("\\'", "", filename$Paper.ID)
filename$Paper.ID <- gsub(" ", "", filename$Paper.ID)

#find number of unique paper IDs per cell and delete row if equal to 1
for (each in c(length(LEEs):1)){
  if (length(unique(unlist(strsplit(filename$Paper.ID[each], ",")))) == 1){
    filename <- filename[-c(each),]
  }
}

