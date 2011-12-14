#/bin/sh
#
#  This sed script replaces 3 character exponents 
#  starting with 0 with 2 characters
#     e-0xx -> e-xx
#     e0xx -> exx
#     E-0xx -> E-xx
#     E0xx -> Exx
#         where
#           x is a digit
#
#   It takes one argument, the file to be operated on.
#   And, it writes to standard out. It may be used to do a 
#   replacement in place.
#
cp $1 .exp.txt
cat .exp.txt | sed 's/\([eE]-\)\(0\)\([0-9][0-9]\)/\1\3/g' | \
               sed 's/\([eE]\)\(0\)\([0-9][0-9]\)/\1\3/g'  | \
               sed 's/\([eE]+\)\(0\)\([0-9][0-9]\)/\1\3/g' 
rm .exp.txt
