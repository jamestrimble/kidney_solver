# Changes each edge weight in a .input or .ndds file to a
# random number in the range [0, 1]
BEGIN { OFS = "\t"; }
NR == 1 { print; }
$1 == -1 { print; exit; }
NR > 1 { print $1, $2, rand(); }
