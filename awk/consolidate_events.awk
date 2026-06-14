#
# count all events in a cevents.dat file that happen 
# at same time as one event
#
# $Id: consolidate_events.awk,v 1.1 2003/01/17 01:12:02 becker Exp $
#
#
BEGIN{
  old_time=-1;
}
{
  time=$1;
  moment=$3;

  if(old_time == -1){
    old_time=time;
    sum_moment=moment;
  }else{
    if(time == old_time){
      sum_moment += moment;
    }else{
      print(old_time,sum_moment);
      sum_moment=moment;
      old_time=time;
    }
  }
}
