#
# compute mean values for all column
#
BEGIN{
}
{
  if((substr($1,1,1)!="#")){
    if(NF > nfmax)
      nfmax=NF;
    for(i=1;i <= NF;i++){
      if(tolower($i) != "nan"){
	sum[i] += $i;
	n[i]++;
      }
    }
  }
}
END{
  for(i=1;i <= nfmax;i++){
    if(n[i]==0)
      printf("NaN ");
    else
      printf("%22.16e ",sum[i]/n[i]);
  }
  printf("\n");
}
