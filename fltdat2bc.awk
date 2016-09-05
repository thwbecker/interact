#
# convert a flt.dat type output of fault patch to a bc.in
# type input file
#
BEGIN{
  useslip=0;    # all slips, if non-zero
  use2stress=1; # strike and normal stress only
    
}
{
  if(substr($1,1,1)!="#"){
    if(useslip){
      if($6!=0) # strike
	print($12,0,$6);
      if($7!=0) # normal
	print($12,2,$7);
      if($8!=0) # dip
	print($12,1,$8);
      
    }
    if(use2stress){
      printf("%i %i %15.6e\n",$12,10,-$9);
      printf("%i %i %15.6e\n",$12,30,-$10);
    }
    

  }
}
