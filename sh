#!/usr/local/bin/bash-5.1.8

declare -A args

args[t]=t
args[g]=top_upperleft
args[p]=p
args[f_cdfpos]=F_pos
args[g_cdfpos]=G_pos
args[f_cdfneg]=F_neg
args[g_cdfneg]=G_neg
args[ff]=f
args[gg]=g
args[fpos]=fpos
args[fneg]=fneg
args[gpos]=gpos
args[gneg]=gneg
args[qfneg]=qfneg
args[qfpos]=qfpos
args[qgpos]=qgpos
args[qgneg]=qgneg

input_name=top_syn
output_name=tmp

full_input_ref=wavz_synthetic
full_input_ref_top=top_syn

full_input_test=wavz_0.3000_0.3000
full_input_test_top=top_upperleft

full_args=""

for key in ${!args[@]}
    do full_args="$full_args$key=${args[$key]}.rsf "
done

echo "full_args = $full_args"

sfwindow n1=1 f1=0 n2=1 f2=10 < $full_input_ref.rsf | sfput n1=1 n2=1 n3=$1 > $full_input_ref_top.rsf
sfwindow n1=1 f1=0 n2=1 f2=10 < $full_input_test.rsf | sfput n1=1 n2=1 n3=$1 > $full_input_test_top.rsf

#./ot.exe g=top_upperleft.rsf t=t.rsf p=p.rsf f_cdfpos=F_pos.rsf f_cdfneg=F_neg.rsf g_cdfpos=G_pos.rsf g_cdfneg=G_neg.rsf ff=f.rsf fneg=fneg.rsf fpos=fpos.rsf gg=g.rsf gpos=gpos.rsf gneg=gneg.rsf qfpos=qfpos.rsf qfneg=qfneg.rsf qgpos=qgpos.rsf qgneg=qgneg.rsf < top_syn.rsf > tmp.rsf

./ot.exe $full_args < $input_name.rsf > $output_name.rsf

(
   create_rsf(){
      sfput n1=$1 n2=1 n3=1 < $2.rsf | sfgraph title="$2" > $2.vpl 
   }

   sfp(){
      sfpen $1.vpl &
   }

   cp(){
       echo "Creating RSF and VPL for ${args[$2]}"
       create_rsf $1 ${args[$2]}
       sfp ${args[$2]}
   }
   
   cp $1 f_cdfpos
   cp $1 f_cdfneg
   cp $1 g_cdfpos
   cp $1 g_cdfneg
   cp $1 ff
   cp $1 fpos
   cp $1 fneg
   cp $1 gg
   cp $1 gpos
   cp $1 gneg
   cp $1 qfneg
   cp $1 qfpos
   cp $1 qgpos
   cp $1 qgneg
)

#sfput n1=$1 n2=1 n3=1 < F_pos.rsf | sfgraph title="F_pos" > F_pos.vpl
#sfput n1=$1 n2=1 n3=1 < F_neg.rsf | sfgraph title="F_neg" > F_neg.vpl
#sfput n1=$1 n2=1 n3=1 < G_pos.rsf | sfgraph title="G_pos" > G_pos.vpl
#sfput n1=$1 n2=1 n3=1 < G_neg.rsf | sfgraph title="G_neg" > G_neg.vpl
#sfput n1=$1 n2=1 n3=1 < f.rsf | sfgraph title="f" > f.vpl
#sfput n1=$1 n2=1 n3=1 < fpos.rsf | sfgraph title="fpos" > fpos.vpl
#sfput n1=$1 n2=1 n3=1 < fneg.rsf | sfgraph title="fneg" > fneg.vpl
#sfput n1=$1 n2=1 n3=1 < g.rsf | sfgraph title="g" > g.vpl
#sfput n1=$1 n2=1 n3=1 < gpos.rsf | sfgraph title="gpos" > gpos.vpl
#sfput n1=$1 n2=1 n3=1 < gneg.rsf | sfgraph title="gneg" > gneg.vpl
#sfput n1=$2 n2=1 n3=1 < qfpos.rsf | sfgraph title="Quantile F_pos" > qfpos.vpl
#sfput n1=$2 n2=1 n3=1 < qfneg.rsf | sfgraph title="Quantile F_neg" > qfneg.vpl
#sfput n1=$2 n2=1 n3=1 < qgpos.rsf | sfgraph title="Quantile G_pos" > qgpos.vpl
#sfput n1=$2 n2=1 n3=1 < qgneg.rsf | sfgraph title="Quantile G_neg" > qgneg.vpl
#
#sfpen F_pos.vpl &
#sfpen G_pos.vpl &
#sfpen F_neg.vpl &
#sfpen G_neg.vpl &
#sfpen f.vpl &
#sfpen fpos.vpl & 
#sfpen fneg.vpl & 
#sfpen g.vpl & 
#sfpen gpos.vpl &
#sfpen gneg.vpl &
#sfpen qfpos.vpl &
#sfpen qfneg.vpl &
#sfpen qgpos.vpl &
#sfpen qgneg.vpl &

A=$(sfdisfil < $output_name.rsf)

echo $A

