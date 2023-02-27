#!/bin/bash

# This is the input: you can use $1 and $2 to read input as cmd line argument
FILE=$1 #'bash_hist_test.dat'
CHANNEL_NUMBER=$2  # They are actually 10: 0 is already a channel
DECIMALS=3

# check the max and the min to define the dimension of the channels:
MAX=`sort -g $FILE | tail -n 1`
MIN=`sort -rg $FILE | tail -n 1`
echo "max=" $MAX
echo "min=" $MIN

# Define the channel width 
CHANNEL_DIM_LONG=`echo "($MAX-$MIN)/($CHANNEL_NUMBER)" | bc -l`
echo "channel dim long = " $CHANNEL_DIM_LONG
#CHANNEL_DIM=`printf '%2.2f' $CHANNEL_DIM_LONG`
CHANNEL_DIM=`echo "{scale=$DECIMALS;print $CHANNEL_DIM_LONG/1;}" | bc -l`
echo "channel dim = " $CHANNEL_DIM

# Probably printf is not the best function in this context because
#+the result could be system dependent.

# Determine the channel for a given number
# Usage: find_channel <number_to_histogram> <width_of_histogram_channel>
function find_channel(){
  NUMBER=$1
  CHANNEL_DIM=$2

  # The channel is found dividing the value for the channel width and 
  #+rounding it.
  RESULT_LONG=`echo $NUMBER/$CHANNEL_DIM | bc -l`
  RESULT=`echo "{scale=0;print $RESULT_LONG/1;}" | bc -l`
  #RESULT=`printf '%.0f' $RESULT_LONG`
  echo $RESULT
}

# Read the file and do the computuation
while IFS='' read -r line || [[ -n "$line" ]]; do

  CHANNEL=`find_channel $line $CHANNEL_DIM`

  [[ -z HIST[$CHANNEL] ]] && HIST[$CHANNEL]=0
  let HIST[$CHANNEL]+=1
done < $FILE

counter=0
for i in ${HIST[*]}; do
  CHANNEL_START=`echo "$CHANNEL_DIM * $counter - .04" | bc -l`
  CHANNEL_END=`echo " $CHANNEL_DIM * $counter + .05" | bc`
  echo "{scale=1;print $CHANNEL_START, \":\", $CHANNEL_END,  \"=>\", $i, \"\n\";}" | bc -l
  #printf '%+2.1f : %2.1f => %i\n' $CHANNEL_START $CHANNEL_END $i
  let counter+=1
done
