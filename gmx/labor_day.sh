
#-----------------------------------------no hoh, 100-200ps, 10ps间隔---------------------------------------------------------------------------------

script_dir="/home/wurp/PycharmProjects/BioUtil/gmx/pygmx"
work_dir="/home/wurp/workspace/antibody/SARS-COV-2"

 anti_list=("1_7KFY" "2_7KFX" "3_7KFV" "4_7KFW" "5_7JVA" "6_7KGK" "7_6LZG" "8_6YZ5" "9_6ZBP" "10_7B27" "11_7BWJ")
 ri_left=(196    196     196     196      1           196       1          196    196    203    195)
 ri_right=(632    634    634     631      227         314       596         322    323    327    636)
 li_left=(1      1        1       1        228          1       597        1      1      1      1)
 li_right=(195   195     195     195      420          195      791        195    195    202    194)

#anti_list=("12_7CH4" "13_7CH5" "14_7E23" "15_7JMO" "16_7K8M" "17_6W41" "18_6YM0" "19_6ZER" "20_7C01" "21_7DEO" "22_7MZF" "23_7DPM")
#ri_left=(   1          1         195       194       1          1        198       198       196       1         197       1       )
#ri_right=(  431        424       403       622       428        443      634       633       628       241       624       444     )
#li_left=(   432        425       1         1         429        444      1         1         1         242       1         445     )
#li_right=(  619        607       194       193       612        638      197       197       195       443       196       640     )

# anti_list=("13_7CH5" "14_7E23" "15_7JMO" "16_7K8M" "18_6YM0" "19_6ZER" "20_7C01" "21_7DEO" "22_7MZF" "23_7DPM")
# ri_left=(   1         195       194       1         198       198       196       1         197       1       )
# ri_right=(  424       403       622       428       634       633       628       241       624       444     )
# li_left=(   425       1         1         429       1         1         1         242       1         445     )
# li_right=(  607       194       193       612       197       197       195       443       196       640     )


cnt=0

for anti in "${anti_list[@]}"
do
    this_dir=$work_dir/$anti/MD_10ns
    cd "$this_dir"
    # mkdir -p most
    # cd most
    # python $script_dir/main.py -tpr "$this_dir"/md_0.tpr -xtc "$this_dir"/md_0.xtc -ri \
    # "${ri_left[$cnt]}" "${ri_right[$cnt]}" -li "${li_left[$cnt]}" "${li_right[$cnt]}" -fm most -rm normal -t 1000 5000 20 1
    mkdir -p 1-5-20
    cd 1-5-20
    python $script_dir/main.py -tpr "$this_dir"/md_0.tpr -xtc "$this_dir"/md_0.xtc -ri \
    "${ri_left[$cnt]}" "${ri_right[$cnt]}" -li "${li_left[$cnt]}" "${li_right[$cnt]}" -fm most -rm normal -t 1000 5000 20 1
    cnt=$cnt+1
done