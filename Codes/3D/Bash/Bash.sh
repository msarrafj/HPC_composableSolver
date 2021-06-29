# Change filename variable and <problem>.py
# nx for this analysis is selected such that dof roughly doubles at each step
# Use this Guidline:
# For CG (TET/HEX): nx = 13 17 21 28 36 45 57
# For DG (TET): nx = 5 6 8 10 13 16 20 
# For DG (HEX): nx = 7 9 11 14 19 23 29 
# For RT0 (TET): nx = 8 10 13 17 22 28 35
# For RT0 (HEX): nx = 13 17 22 28 37 46 58
 
declare -a tolerance=(1e-5 1e-7)
declare -a solver_type=("fields" "scales")
declare -a element=("HEX")
declare -a nx=(7 9 11 14 19 23 29)

mkdir -p ./Output/Convergence ./Output/DoE ./Output/Logs ./Output/Static_scaling ./Output/True_scaling

for l in "${tolerance[@]}"
 do
 echo "------------Tolerance is: $l---------------"


 for k in "${solver_type[@]}"
 do
 echo "Solver_type is: $k"


  for j in "${element[@]}"
  do
  echo "element is: $j"

	filename="DG_${k}split_${j}_tol${l}"
	echo "**********Solve for $filename***********"
	#mv ./$filename.txt ./old_$filename.txt
	for i in "${nx[@]}"  
	do 
	echo " The nx is: ${i}"       
	echo "===================== The nx is: ${i}  =====================" >> ./$filename.txt
	mpiexec -n 1 python ../3D_composable_solvers_DG.py ${l} 0 $k hypre ${i} $j  >> ./$filename.txt
	done
	
	# CONVERGENCE RESULTS
	awk '$1=="DoS" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="DoA_p1" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/Convergence/p1_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat

	awk '$1=="DoS" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="DoA_p2" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/Convergence/p2_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat

	awk '$1=="DoS" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="DoA_v1" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/Convergence/v1_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat
	
	awk '$1=="DoS" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="DoA_v2" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/Convergence/v2_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat
	
	# STATIC_SCALING RESULTS
	awk '$1=="TotaltimeLOG" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="DofSecLOG" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/Static_scaling/$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat

	# TRUE_SCALING
	awk '$1=="TotaltimeLOG" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="truedof_p1" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/True_scaling/p1_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat

	awk '$1=="TotaltimeLOG" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="truedof_p2" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/True_scaling/p2_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat

	awk '$1=="TotaltimeLOG" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="truedof_v1" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/True_scaling/v1_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat
	
	awk '$1=="TotaltimeLOG" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="truedof_v2" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/True_scaling/v2_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat

	# DoE RESULTS
	awk '$1=="TotaltimeLOG" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="DoE_p1" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/DoE/p1_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat

	awk '$1=="TotaltimeLOG" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="DoE_p2" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/DoE/p2_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat

	awk '$1=="TotaltimeLOG" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="DoE_v1" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/DoE/v1_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat
	
	awk '$1=="TotaltimeLOG" {print $3}' ./$filename.txt > ./results1_$filename.dat
	awk '$1=="DoE_v2" {print $3}' ./$filename.txt > ./results2_$filename.dat
	paste ./results1_$filename.dat ./results2_$filename.dat > ./Output/DoE/v2_$filename.dat
	rm ./results1_$filename.dat ./results2_$filename.dat

	# MOVE filename.txt to Output/Logs dir
	mv ./$filename.txt ./Output/Logs/.	
  done
 done
done
