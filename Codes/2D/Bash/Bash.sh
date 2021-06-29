# Change filename variable and <problem>.py 
declare -a tolerance=(1e-5 1e-6 1e-7 1e-8)
declare -a solver_type=("fields" "scales")
declare -a element=("T3" "Q4")
declare -a nx=(5 10 20 40 80 160)

mkdir -p ./Output/Convergence ./Output/DoE ./Output/Logs ./Output/Static_scaling ./Output/True_scaling

for l in "${tolerance[@]}"
 do
 echo "------------Tolerance is: $l---------------"


 for k in "${solver_type[@]}"
 do
 echo "Solver_type is: $k"


  for j in "${element[@]}"
  do
  echo "amg_type is: $j"

	filename="RT0_${k}split_${j}_tol${l}"
	echo "**********Solve for $filename***********"
	#mv ./$filename.txt ./old_$filename.txt
	for i in "${nx[@]}"  
	do 
	echo " The nx is: ${i}"       
	echo "===================== The nx is: ${i}  =====================" >> ./$filename.txt
	mpiexec -n 1 python ../2D_composable_solver_RT0.py ${l} 0 $k hypre ${i} $j  >> ./$filename.txt
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
