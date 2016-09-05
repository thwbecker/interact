#
# read cevents.dat and print aki style focal mechanism
#
BEGIN{

  pif=57.295779513082320876798154814105;
  
  d=50;offa=0.;

  doffa=20/pif;

  n=0;
}
{

  n++;
  time=$1;
  mom=$3;

  us=$5;ud=$6;			# strike, dip slip

  u[1]=$8;u[2]=$9;u[3]=$10;	# cartesian slip
  
  x=$11;y=$12;depth=-$13;

  mag=log(mom)-7;

  # offset for plotting

  x0 = x + d*cos(offa);
  y0 = y + d*sin(offa);
  offa+=doffa;

  title=sprintf("%.2f",time);
  
  strike=atan2(u[1],u[2])*pif;
  dip = 90-atan2(u[3],sqrt(u[1]**2+u[2]**2))*pif;
  if(dip>90){
    dip = 180-dip;
    strike += 180;
  }
  if(dip<0){
    dip = -dip;
    strike+=180;
  }
  if(strike<0)
    strike+=360;
  if(strike>360)
    strike-=360;

  rake= atan2(-ud,us)*pif;

#     (a) Focal mechanism in AKI & RICHARD's convention:
#          X, Y, depth, strike, dip, rake, mag, newX, newY, event_title
  print(x,y,depth,strike,dip,rake,mag,x0,y0,title)
#  print(x,y,x0,y0)
}