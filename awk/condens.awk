#
# "condens" several events in an ascii event file to 
# one group event, sum individual moments
#
# $Id: condens.awk,v 1.1 2003/01/17 01:11:32 becker Exp $
#
BEGIN{
  tm=0;
}
{
  if(NR==1){
    t=$1;
    oldgroup=$3;
    moment=$4;
  }else{
    if($1 == t){
      if($3==oldgroup)
	moment += $4;
      else{	
	tm += moment;
	print(t,oldgroup,moment,tm);
	t=$1;
	moment=$4;
	oldgroup=$3;
      }
    }else{
      tm += moment;
      print(t,oldgroup,moment,tm);
      t=$1;
      moment=$4;
      oldgroup=$3;
    }
  }
}
END{
  tm += moment;
  print(t,oldgroup,moment,tm);

}
