
# /media/xin/WinData/ACS/github/BioUtil/gmx/Loop_MD.sh


# anti_list=("1_7KFY" "2_7KFX" "3_7KFV" "4_7KFW" "5_7JVA" "6_7KGK" "7_6LZG" "8_6YZ5" "9_6ZBP" "10_7B27" "11_7BWJ")
anti_list=("12_7CH4" "13_7CH5" "14_7E23" "15_7JMO" "16_7K8M" "17_6W41" "18_6YM0" "19_6ZER" "20_7C01" "21_7DEO" "22_7MZF" "23_7DPM")
# anti_list=("12_7CH4")


for anti in ${anti_list[@]}
do
    # anti_path=/home/wurp/workspace/antibody/SARS-COV-2/$anti
    anti_path=/media/xin/Raid0/ACS/gmx/interaction/$anti
    /media/xin/WinData/ACS/github/BioUtil/gmx/MD.sh $anti_path

       
done


#  /home/wurp/workspace/antibody/script/MD_batchRun.sh