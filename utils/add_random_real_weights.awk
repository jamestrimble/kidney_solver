# Changes each edge weight in a .input or .ndds file to a
# random number in the range [0, 1)

# Optionally, use -v seed=$RANDOM from the command line
# to set the random-number generator's seed randomly

BEGIN {
    srand(seed ? seed : 1); 
    OFS = "\t";
}

# If the line represents and edge, replace weight with a random number 
NF == 3 && $1 != -1 { print $1, $2, rand(); next }

# Otherwise, just print the line without changes
{ print }
