#
# prints the norm of entries in a line
#
BEGIN{
  
}
{
  if((substr($1,1,1)!="#")){
    sum=0.0;
    for(i=1;i<=NF;i++){
      if(tolower($i)!="nan")
	sum += ($i)*($i);
    }
    print(sqrt(sum));
  }
}

