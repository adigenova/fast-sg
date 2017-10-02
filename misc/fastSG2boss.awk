
BEGIN{OFS="\t"; i=1;}
{
 if($0 ~/^@/){
	print $0 > name".fwd.sam";
	print $0 > name".rev.sam";
 }else{ 
       $1=$1"/1";
       print $0 > name".fwd.sam";
       getline;
       $1=$1"/2";
       print $0 > name".rev.sam"; 
       i++;
 }
}

