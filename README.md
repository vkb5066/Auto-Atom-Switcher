Auto atom switcher - updated  8/12/19

Input:	1. POSCAR-style file, located in directory defined in the PoscarInfo.h header file, on lines 15 and 16
	   	2. bond data info file, located in directory defined in the PoscarInfo.h header file, on lines 15 and 17

Output: 	1. User-defined number of edited versions of the input POSCAR-style file
			-The C++ code defaults to placing the new files to the input file directory, but the provided shell script moves these files 
			 into a sub-directory named 'poscars'
		2. (optional) A .csv file with information on the generated files
		3. (optional) A log file with information of every defect found per iteration, for each generated file
			-The C++ code defaults to placing these log files to the input file directory, but the provided shell script moves these files
			into a sub-directory named 'poscars/defectInfo'

To Use:
		1. Edit PoscarInfo.h's lines 15 - 17 to the desired file paths
		2. Compile the .cpp file with gcc
			-The shell script assumes that the compiled file is called 'drive', so you may want to change the name from 'a.out' to 'drive'
		3. Run the shell script, making sure that there is no directory named 'poscars' in the current directory (else it will be deleted without 			warning)
		4. Follow instruction printed to the terminal