#/bin/sh


if [ "$1" == "1" ];
then
    obj="onemax"
fi
if [ "$1" == "2" ];
then
    obj="needle"
fi
if [ "$1" == "3" ];
then
    obj="dectrap"
fi
if [ "$1" == "4" ];
then
    obj="dec2trap"
fi
if [ "$1" == "5" ];
then
    obj="2trap"
fi
if [ "$1" == "6" ];
then
    obj="royal_road"
fi
if [ "$1" == "7" ];
then
    obj="longpath"
fi



mkdir data
rm -rf "data/$obj.$2pop"
mkdir "data/$obj.$2pop"
    

#==========---------------
echo "rechenberg_unitation - $obj - $2 pop member"
for i in `seq 1 25`;
do
   adaptive/rechenberg_unitation $1 $2 > junk
   mv bog.dat data/$obj.$2pop/rechenberg_unitation_$obj_$2pop_$i.dat
   echo $i
done 


