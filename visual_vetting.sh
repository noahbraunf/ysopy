# Purpose: neat sh program to display images in same window, close them, and allow user to classify them

sudo apt-get install fim
for i $(ls)
  do 
  fim $i
  read varname
  echo $i, $varname
  done
 
