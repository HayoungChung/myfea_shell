
for T in 0.001 0.01 0.1 1
do
mkdir $T
cp optim_base $T
cd $T
./optim_base $T
cd ..
done

