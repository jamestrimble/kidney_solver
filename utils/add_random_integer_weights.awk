# Set the following variables on the command line:
# lower, upper

# Changes each edge weight in a .input or .ndds file to a
# random integer number in the range [lower, upper)

# Optionally, use -v seed=$RANDOM from the command line
# to set the random-number generator's seed randomly

BEGIN {
    if (length(upper)==0 || length(lower)==0) {
        print "Please specify lower and upper limits; e.g." > "/dev/stderr";
        print "awk -f <scriptname> -v lower=5 -v upper=10" > "/dev/stderr";
        exit 1;
    }
    srand(seed ? seed : 1); 
    OFS = "\t"
}

function randint() {
    return int(lower + (upper-lower) * rand());
}

# If the line represents and edge, replace weight with a random number 
NF == 3 && $1 != -1 { print $1, $2, randint(); next }

# Otherwise, just print the line without changes
{ print }
