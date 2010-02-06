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
    obj="royal road"
fi
if [ "$1" == "7" ];
then
    obj="longpath"
fi



mkdir data
rm -rf "data/$obj.$2pop"
mkdir "data/$obj.$2pop"
    

#==========---------------

echo "static_unitation - $obj - $2 pop member"
for i in `seq 1 25`;
do
   dynamic/static_unitation $1 $2 > junk
   mv bog.dat data/$obj.$2pop/static_unitation_$obj_$2pop_$i.dat
   echo $i
done 


echo "droste_unitation - $obj - $2 pop member"
for i in `seq 1 25`;
do
   dynamic/droste_unitation $1 $2 > junk
   mv bog.dat data/$obj.$2pop/droste_unitation_$obj_$2pop_$i.dat
   echo $i
done 

echo "back_unitation - $obj - $2 pop member"
for i in `seq 1 25`;
do
   dynamic/back_unitation $1 $2 > junk
   mv bog.dat data/$obj.$2pop/back_unitation_$obj_$2pop_$i.dat
   echo $i
done 

#================-------------

echo "constantgain_unitation - $obj - $2 pop member"
for i in `seq 1 25`;
do
   adaptive/constantgain_unitation $1 $2 > junk
   mv bog.dat data/$obj.$2pop/constantgain_unitation_$obj_$2pop_$i.dat
   echo $i
done 

echo "decliningadapt_unitation - $obj - $2 pop member"
for i in `seq 1 25`;
do
   adaptive/decliningadapt_unitation $1 $2 > junk
   mv bog.dat data/$obj.$2pop/decliningadapt_unitation_$obj_$2pop_$i.dat
   echo $i
done 

echo "rechenberg_unitation - $obj - $2 pop member"
for i in `seq 1 25`;
do
   adaptive/rechenberg_unitation $1 $2 > junk
   mv bog.dat data/$obj.$2pop/rechenberg_unitation_$obj_$2pop_$i.dat
   echo $i
done 


#================-------------


#echo "lee_unitation - $obj - $2 pop member"
#for i in `seq 1 25`;
#do
#   fuzzy/lee/lee_unitation $1 $2 > junk
#   mv bog.dat data/$obj.$2pop/lee_unitation_$obj_$2pop_$i.dat
#   echo $i
#done 

echo "shi_unitation - $obj - $2 pop member"
for i in `seq 1 25`;
do
   fuzzy/shi3/shi_unitation $1 $2 > junk
   mv bog.dat data/$obj.$2pop/shi_unitation_$obj_$2pop_$i.dat
   echo $i
done 

echo "Done"

