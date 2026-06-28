/hmat_storage backend/ && b=="" {
  for(i=1;i<=NF;i++){
    if($i=="backend")        b=$(i+1)
    if($i=="m")              N=$(i+1)
    if($i=="stored_scalars") s=$(i+1)
    if($i=="mbytes")         mb=$(i+1)
  }
}
/compression ratio:/ && ratio=="" {ratio=$NF}                 # htool: stored is NA, recover from ratio
/interaction matrix build wall/ && build=="" {build=$(NF-1)}  # assembly cost [s]
/^Time \(sec\):/ {ttot=$3}                                    # total wall [s]
$1=="MatMult" {mmt=$4; mmn=$2}                                # matvec total time, count
END{
  if(s=="NA" && ratio>0){ s=int(N*N/ratio); mb=s*8/1048576 }
  printf "%s %s %s %.1f %s %s %.4f\n", b, N, s, mb, build, ttot, (mmn>0?1000*mmt/mmn:0)
}
