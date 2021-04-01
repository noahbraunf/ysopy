# Purpose: neat sh program to display images in same window, close them, and allow user to classify them

sudo apt-get install fim
for i in 1 2 3 4 5 6 6 #replace these numbers with an array of ids
  do 
  fim #image path
  read varname
  echo $i, $varname
  done
 
