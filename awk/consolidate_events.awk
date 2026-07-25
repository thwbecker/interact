#
# count all events in a cevents.dat file that happen 
# at same time as one event
#
# part of interact, (C) Thorsten Becker; see README.md and COPYRIGHT
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
