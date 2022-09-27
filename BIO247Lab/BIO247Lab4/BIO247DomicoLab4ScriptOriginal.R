#Script for Removing Junk LEEs: ORIGINAL
#Whitney Domico BIO247 Lab 4

#test comment for git
#second test comment

#uploading and defining dataframe
Lab4Data <- read.csv("C:\\Users\\Whitn\\OneDrive\\Desktop\\BIO247\\BIO247Lab\\BIO247Lab4\\MachineRead_output.csv", header=TRUE)

#create vector of paper IDs
LEEs <- Lab4Data$Paper.ID

#remove "junk" characters
Lab4Data$Paper.ID <- gsub("\\]", "", Lab4Data$Paper.ID)
Lab4Data$Paper.ID <- gsub("\\[", "", Lab4Data$Paper.ID)
Lab4Data$Paper.ID <- gsub("\\'", "", Lab4Data$Paper.ID)
Lab4Data$Paper.ID <- gsub(" ", "", Lab4Data$Paper.ID)

#find number of unique paper IDs per cell and delete row if equal to 1
for (each in c(length(LEEs):1)){
  if (length(unique(unlist(strsplit(Lab4Data$Paper.ID[each], ",")))) == 1){
    Lab4Data <- Lab4Data[-c(each),]
  }
}
