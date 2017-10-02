head -n 10000 $1 > $1.1000.txt
diff -s $1.1000.txt example/$1.1000.txt
rm -f $1.1000.txt
