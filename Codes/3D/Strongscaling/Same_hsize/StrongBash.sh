#!/bin/bash
# When running remember that nx should be set to just one number
# file.py changes
# change filename

declare -a MPI=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
declare -a solver_type=("fields" "scales")
declare -a element=("TET" "HEX")
declare -a nx=(16)

TimeToSol=()

for k in "${solver_type[@]}"
  do
  echo "Solver_type is: $k"

  for j in "${element[@]}"
    do
    echo "element is: $j"
	mpiexec -n 1 python ../../3D_composable_solvers_RT0.py 1e-7 1 $k hypre 16 $j
	echo "pre-run....."

     for l in "${MPI[@]}"
      do
      echo "# MPI proc is: ${l}"
      mkdir -p ./MPI_${l}/Convergence ./MPI_${l}/DoE ./MPI_${l}/Logs ./MPI_${l}/Static_scaling ./MPI_${l}/True_scaling ./MPI_${l}/KSP_iterations ./Strong_results ./MPI_efficiency_results
		filename="RT0_${k}split_${j}"
		echo "***********${filename}**********"
		#If the logfile already exists, them move it to old_logfile
		#if [ -f ./MPI_${l}/Logs/$filename.txt ]; then
		#	echo "Old $filename.txt backed up"
		#	mv ./MPI_${l}/Logs/$filename.txt ./MPI_${l}/Logs/old_$filename.txt
		#fi
		# SOLVING THE PROBLEM
		echo " <<${l} MPI processors used in this problem>>" >> ./MPI_${l}/$filename.txt
		for i in "${nx[@]}"  
		do 
		echo "===================== The nx is: ${i}  =====================" >> ./MPI_${l}/$filename.txt
		mpiexec -n ${l} python ../../3D_composable_solvers_RT0.py 1e-7 1 $k hypre ${i} $j  >> ./MPI_${l}/$filename.txt
		done
		# Comment this part if <show_monitor> is 0 (monitor is off)
		awk '$1=="Linear" {print $NF}' ./MPI_${l}/$filename.txt > ./MPI_${l}/KSP_iterations/$filename.dat
		
		# CONVERGENCE RESULTS
		awk '$1=="DoS" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="DoA_p1" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/Convergence/p1_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat

		awk '$1=="DoS" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="DoA_p2" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/Convergence/p2_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat

		awk '$1=="DoS" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="DoA_v1" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/Convergence/v1_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat
		
		awk '$1=="DoS" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="DoA_v2" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/Convergence/v2_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat
		# STATIC_SCALING RESULTS
		awk '$1=="TotaltimeLOG" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="DofSecLOG" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/Static_scaling/$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat
		
		# TRUE_SCALING
		awk '$1=="TotaltimeLOG" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="truedof_p1" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/True_scaling/p1_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat
	
		awk '$1=="TotaltimeLOG" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="truedof_p2" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/True_scaling/p2_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat

		awk '$1=="TotaltimeLOG" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="truedof_v1" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/True_scaling/v1_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat
		
		awk '$1=="TotaltimeLOG" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="truedof_v2" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/True_scaling/v2_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat
	
		# DoE RESULTS
		awk '$1=="TotaltimeLOG" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="DoE_p1" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/DoE/p1_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat

		awk '$1=="TotaltimeLOG" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="DoE_p2" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/DoE/p2_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat

		awk '$1=="TotaltimeLOG" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="DoE_v1" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/DoE/v1_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat
		
		awk '$1=="TotaltimeLOG" {print $3}' ./MPI_${l}/$filename.txt > ./results1_$filename.dat
		awk '$1=="DoE_v2" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat > ./MPI_${l}/DoE/v2_$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat
		
		#STRONG-SCALING OPERATION
		zaman=$(awk '$1=="Totaltime" {print $3}' ./MPI_${l}/$filename.txt)
		#echo "Totaltime is: ${zaman}"
		TimeToSol+=("${zaman}")
		echo "Time to solution is ${TimeToSol[@]}"
		#TimeToSol1=${TimeToSol[0]}
		#echo "$TimeToSol1"
		efficiency=$(awk "BEGIN {print (100 * ${TimeToSol[0]} )/($l * ${zaman})}") 
		echo " Efficiency is: $efficiency"
	
		# lists <MPI>	<efficiency>
		echo "$l" > ./results1_$filename.dat
		echo "$efficiency" > ./results2_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat >> ./MPI_efficiency_results/$filename.dat
		rm ./results1_$filename.dat ./results2_$filename.dat

		# list <MPI> <Assemblytime> <Solvetime> <Totaltime> <ksp_iters> <efficiency>
		echo "$l" > ./results1_$filename.dat
		awk '$1=="Assemblytime" {print $3}' ./MPI_${l}/$filename.txt > ./results2_$filename.dat
		awk '$1=="Solvetime" {print $3}' ./MPI_${l}/$filename.txt > ./results3_$filename.dat
		awk '$1=="Totaltime" {print $3}' ./MPI_${l}/$filename.txt > ./results4_$filename.dat
		awk '$1=="Linear" {print $NF}' ./MPI_${l}/$filename.txt > ./results5_$filename.dat
		echo "$efficiency" > ./results6_$filename.dat
		paste ./results1_$filename.dat ./results2_$filename.dat ./results3_$filename.dat ./results4_$filename.dat ./results5_$filename.dat ./results6_$filename.dat >> ./Strong_results/$filename.dat 
		rm ./results{1..6}_$filename.dat
		 

		# MOVE DATA TO Logs
		mv ./MPI_${l}/$filename.txt ./MPI_${l}/Logs
        done	
		#RESET TIME ARRAY
		TimeToSol=()				
    done
done
