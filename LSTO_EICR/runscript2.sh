
for T in 0.001 0.01 0.1 1
do
mkdir adj_$T
cp optim_base_adj adj_$T
cd adj_$T
./optim_base_adj $T
cd ..
done

