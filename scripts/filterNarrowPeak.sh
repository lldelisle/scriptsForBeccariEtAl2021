if [ $1 = "-h" ] || [ $1 = "--help" ]; then
  echo "Usage: bash filterNarrowPeak.sh input.narrowPeak filtered.narrowPeak"
  exit 0
fi

cat $1 | awk '
BEGIN{
  best=""  # Initialize the best value to ""
}
!($4~/[a-z]$/){  # If the peak name does not end with a letter, no problem
  if(best!=""){  # If there was a best value stored
    print best  # Print it
    best=""  # Initialize the best value to ""
  }
  print  # Print the peak with no letter
}
$4~/[a-z]$/{ # If the peak name ends with a letter
  if($4~/a$/){ # If it ends with a 'a', this is a new peak
    if(best!=""){  # If there was a best value stored
      print best  # Print it
    }
    best=$0  # Store the current peak (the one with 'a') in best
    bestScore=$5  # And its score in bestScore
  }else{  # If the peak name does not start with 'a'
    if($5>bestScore){  # Check if it is better than the one stored
      bestScore=$5  # If it is the case, change the bestScore
      best=$0  # And store the best
    }
  }
}' > $2
